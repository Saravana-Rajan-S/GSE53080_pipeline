#!/usr/bin/env python3
"""
GSE53080 — 28 WARNING Files: 5nt vs 9nt Diagnostic v4
======================================================
Repo: src/preprocessing/verify_28_warning_files.py

Purpose:
  Re-scan the 28 WARNING files from barcode verification to determine
  whether the low 9nt adapter rate is caused by:
    A) TCGT ERROR     — barcode is present but the TCGT suffix is degraded
    B) FALSE POSITIVE — 5nt barcode appears twice in same read (once real,
                        once by random chance inside the RNA sequence)

Methodology:
  1. Loads barcode_verification_report.csv to get the pre-assigned barcode
     for each WARNING file (no re-scanning to find dominant adapter).

  2. SINGLE-PASS scan through every read — for each read:
       a) Search for full 9nt signature (exact, 0 mismatch, sliding).
          If found -> 9nt match. Record position (1-based).
          Since 5nt is the prefix of 9nt, 5nt is implicitly present here.
          Use continue — skip 5nt search for this read (no double-counting).

       b) If 9nt NOT found -> count how many times 5nt barcode appears
          in this read using seq.count(barcode_5nt):
            count == 0 -> no match at all
            count == 1 -> barcode found once -> real barcode, TCGT degraded
            count >= 2 -> barcode found 2+ times -> one is real adapter,
                         the other(s) are random occurrences inside RNA.
                         This read contributes to the double_hit_reads count.
          Record the position of the FIRST 5nt occurrence (1-based).

  3. CHECKPOINT: fresh 9nt% must match previous Adapter_Percentage within 1%.
     Confirms this script uses identical scan logic to verify_barcodes.py.

  4. CONCLUSION logic:
       double_hit_pct >= 10% -> FALSE POSITIVE (partial)
         Significant false positive contamination. 9nt is reliable measure.
       double_hit_pct <  10% -> TCGT ERROR confirmed
         Nearly all 5nt matches are genuine adapters with broken TCGT.

  5. TRIMMING: Always trim at the 5nt barcode position regardless of conclusion.
     The diagnostic is for biological understanding only.

Read structure:
  [RNA insert, variable length] + [5nt barcode] + TCGTATGCCGTCTTCTGCTTG
  Trim point: right before the 5nt barcode, wherever it occurs in the read.

Outputs:
  docs/warning_files_5nt_vs_9nt_diagnostic.csv

Prerequisites:
  docs/barcode_verification_report.csv (run verify_barcodes.py first)

Usage:
    cd ~/Documents/GSE53080_pipeline
    python3 src/preprocessing/verify_28_warning_files.py

    # With custom base path:
    BASE_PATH=/path/to/data python3 src/preprocessing/verify_28_warning_files.py
"""

import os
import sys
import gzip
import csv
import time
from pathlib import Path
from datetime import datetime, timedelta

# ============================================================
# CONFIGURATION
# ============================================================
SCRIPT_DIR   = Path(__file__).resolve().parent.parent.parent  # repo root
BASE_PATH    = os.environ.get("BASE_PATH", "/home/saravana/Downloads/GSE53080")
REPORT_PATH  = SCRIPT_DIR / "docs" / "barcode_verification_report.csv"
OUTPUT_PATH  = SCRIPT_DIR / "docs" / "warning_files_5nt_vs_9nt_diagnostic.csv"

# 20 adapters: 5nt barcode + 9nt signature + full 26nt adapter
ADAPTERS = {
    "Ad01": {"barcode": "TCACT", "signature": "TCACTTCGT", "full": "TCACTTCGTATGCCGTCTTCTGCTTG"},
    "Ad02": {"barcode": "TCATC", "signature": "TCATCTCGT", "full": "TCATCTCGTATGCCGTCTTCTGCTTG"},
    "Ad03": {"barcode": "TCCAC", "signature": "TCCACTCGT", "full": "TCCACTCGTATGCCGTCTTCTGCTTG"},
    "Ad04": {"barcode": "TCCGT", "signature": "TCCGTTCGT", "full": "TCCGTTCGTATGCCGTCTTCTGCTTG"},
    "Ad05": {"barcode": "TCCTA", "signature": "TCCTATCGT", "full": "TCCTATCGTATGCCGTCTTCTGCTTG"},
    "Ad06": {"barcode": "TCGAT", "signature": "TCGATTCGT", "full": "TCGATTCGTATGCCGTCTTCTGCTTG"},
    "Ad07": {"barcode": "TCGCG", "signature": "TCGCGTCGT", "full": "TCGCGTCGTATGCCGTCTTCTGCTTG"},
    "Ad08": {"barcode": "TCTAG", "signature": "TCTAGTCGT", "full": "TCTAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad09": {"barcode": "TCTCC", "signature": "TCTCCTCGT", "full": "TCTCCTCGTATGCCGTCTTCTGCTTG"},
    "Ad10": {"barcode": "TCTGA", "signature": "TCTGATCGT", "full": "TCTGATCGTATGCCGTCTTCTGCTTG"},
    "Ad11": {"barcode": "TTAAG", "signature": "TTAAGTCGT", "full": "TTAAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad12": {"barcode": "TAACG", "signature": "TAACGTCGT", "full": "TAACGTCGTATGCCGTCTTCTGCTTG"},
    "Ad13": {"barcode": "TAATA", "signature": "TAATATCGT", "full": "TAATATCGTATGCCGTCTTCTGCTTG"},
    "Ad14": {"barcode": "TAGAG", "signature": "TAGAGTCGT", "full": "TAGAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad15": {"barcode": "TAGGA", "signature": "TAGGATCGT", "full": "TAGGATCGTATGCCGTCTTCTGCTTG"},
    "Ad16": {"barcode": "TATCA", "signature": "TATCATCGT", "full": "TATCATCGTATGCCGTCTTCTGCTTG"},
    "Ad17": {"barcode": "TGATG", "signature": "TGATGTCGT", "full": "TGATGTCGTATGCCGTCTTCTGCTTG"},
    "Ad18": {"barcode": "TGTGT", "signature": "TGTGTTCGT", "full": "TGTGTTCGTATGCCGTCTTCTGCTTG"},
    "Ad19": {"barcode": "TTACA", "signature": "TTACATCGT", "full": "TTACATCGTATGCCGTCTTCTGCTTG"},
    "Ad20": {"barcode": "TTGGT", "signature": "TTGGTTCGT", "full": "TTGGTTCGTATGCCGTCTTCTGCTTG"},
}

ID_TO_SIG     = {k: v["signature"] for k, v in ADAPTERS.items()}
ID_TO_BARCODE = {k: v["barcode"]   for k, v in ADAPTERS.items()}
ID_TO_FULL    = {k: v["full"]      for k, v in ADAPTERS.items()}


# ============================================================
# FUNCTIONS
# ============================================================

def load_warning_files_from_report(report_path):
    """
    Read barcode_verification_report.csv and extract WARNING rows.
    Returns list of dicts: srr_id, folder, assigned_adapter_id,
    assigned_barcode, previous_pct.
    """
    warning_files = []
    with open(report_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("Status", "").strip() == "WARNING":
                try:
                    prev_pct = float(row["Adapter_Percentage"])
                except (ValueError, KeyError):
                    prev_pct = 0.0
                warning_files.append({
                    "srr_id":              row["SRR_ID"],
                    "folder":              row["Folder"],
                    "assigned_adapter_id": row.get("Assigned_Adapter_ID", "").strip(),
                    "assigned_barcode":    row.get("Assigned_Barcode", "").strip(),
                    "previous_pct":        prev_pct,
                })
    return warning_files


def find_fastq(base_path, folder_rel, srr_id):
    """Find the FASTQ file for a given SRR. Tries multiple extensions."""
    folder = Path(base_path) / folder_rel
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        candidate = folder / f"{srr_id}{ext}"
        if candidate.exists():
            return candidate
    return None


def compute_median(values):
    """Compute median. Returns 'N/A' if list is empty."""
    if not values:
        return "N/A"
    sorted_vals = sorted(values)
    n   = len(sorted_vals)
    mid = n // 2
    if n % 2 == 1:
        return sorted_vals[mid]
    else:
        return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2


def scan_file(filepath, barcode_5nt, sig_9nt):
    """
    Single-pass sliding scan through the entire FASTQ file.

    For each read:

      STEP 1 — Search for 9nt signature (exact, 0 mismatch).
        seq.find(sig_9nt) slides through whole read — no fixed position.
        If found:
          - Count it, record 1-based position
          - 'continue' to next read (5nt is prefix of 9nt, already counted)

      STEP 2 — 9nt not found -> use seq.count(barcode_5nt):

        WHY count() not find():
          find() only tells you IF the barcode exists and WHERE (first hit).
          count() tells you HOW MANY TIMES it exists in the read.
          A 5nt barcode (4^5 = 1024 combos) can appear by random chance
          inside the RNA itself. count() exposes this.

        count == 0: no barcode found, skip
        count == 1: barcode found exactly once -> real adapter, TCGT broken
                    Record position of first (only) hit.
        count >= 2: barcode found 2+ times -> one is the real adapter,
                    the extra(s) are random hits inside the RNA sequence.
                    Still counts as a real adapter read (match_5nt_only++),
                    but also flagged as double_hit (false positive signal).
                    Record position of first hit.

    Returns:
      total_reads, match_9nt, pos_9nt_list,
      match_5nt_only, pos_5nt_only_list,
      double_hit_reads, read_lengths
    """
    total_reads       = 0
    match_9nt         = 0
    pos_9nt_list      = []
    match_5nt_only    = 0
    pos_5nt_only_list = []
    double_hit_reads  = 0
    read_lengths      = []

    opener = gzip.open if str(filepath).endswith(".gz") else open

    with opener(filepath, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()

            if not seq or not plus or not qual:
                print(f"    WARNING: Truncated FASTQ at read {total_reads + 1}. Stopping.")
                break

            total_reads += 1
            read_lengths.append(len(seq))

            # STEP 1: Search for full 9nt signature
            pos9 = seq.find(sig_9nt)
            if pos9 != -1:
                match_9nt += 1
                pos_9nt_list.append(pos9 + 1)  # 1-based
                continue                        # skip 5nt search

            # STEP 2: 9nt not found -> count 5nt occurrences
            count_5nt = seq.count(barcode_5nt)

            if count_5nt == 0:
                pass  # no barcode, no match

            elif count_5nt == 1:
                # One hit: real adapter with degraded TCGT
                match_5nt_only += 1
                pos5 = seq.find(barcode_5nt)
                pos_5nt_only_list.append(pos5 + 1)  # 1-based

            else:
                # 2+ hits: real adapter + false positive noise inside RNA
                match_5nt_only += 1          # still a real adapter read
                double_hit_reads += 1        # flag for false positive rate
                pos5 = seq.find(barcode_5nt)
                pos_5nt_only_list.append(pos5 + 1)  # record first hit position

    return (total_reads, match_9nt, pos_9nt_list,
            match_5nt_only, pos_5nt_only_list,
            double_hit_reads, read_lengths)


def build_conclusion(checkpoint_pass, checkpoint_diff,
                     pct_9nt, previous_pct,
                     pct_5nt_total, pct_5nt_only,
                     double_hit_reads, total_reads,
                     pos_9nt_median, pos_5nt_only_median):
    """
    Build Notes and Trimming_Recommendation.

    Decision tree:
      1. Checkpoint fail -> flag immediately, stop.
      2. double_hit_pct >= 10% -> FALSE POSITIVE (partial)
         Significant fraction of 5nt matches are noise from random
         barcode occurrences inside the RNA sequence.
      3. double_hit_pct <  10% -> TCGT ERROR confirmed
         Barcode is present and real in nearly all reads.
         TCGT suffix is what is degraded.

    Threshold of 10%:
      Conservative. A 5nt barcode (1024 combos) in a 51nt read has a
      small but nonzero chance of appearing twice. If >10% of reads
      show double hits, the contamination is systematic, not random chance.

    Trimming is ALWAYS at the 5nt barcode position regardless of conclusion.
    """
    notes_parts  = []
    trimming_rec = "Trim at 5nt barcode position (remove everything from barcode onward)"

    # Checkpoint
    if not checkpoint_pass:
        notes_parts.append(
            f"REPRODUCIBILITY CHECK FAIL: 9nt% changed by {checkpoint_diff:.2f}% "
            f"(previous={previous_pct:.2f}%, fresh={pct_9nt:.2f}%). "
            f"Verify scan logic matches verify_barcodes.py."
        )
        trimming_rec = "Investigate — reproducibility check failed"
        return " | ".join(notes_parts), trimming_rec

    double_hit_pct = (double_hit_reads / total_reads * 100.0) if total_reads > 0 else 0.0
    diff = pct_5nt_total - pct_9nt

    if double_hit_pct >= 10.0:
        notes_parts.append(
            f"FALSE POSITIVE (partial): {double_hit_pct:.2f}% of reads have the 5nt "
            f"barcode appearing 2+ times — one real adapter hit plus random "
            f"occurrence(s) inside the RNA sequence. "
            f"5nt% ({pct_5nt_total:.2f}%) is inflated. "
            f"9nt% ({pct_9nt:.2f}%) is the reliable adapter rate. "
            f"Gap={diff:.2f}% includes both TCGT degradation and false positive noise."
        )
    else:
        notes_parts.append(
            f"TCGT ERROR confirmed: double-hit rate={double_hit_pct:.2f}% (below 10% threshold). "
            f"Nearly all 5nt matches are genuine adapters. "
            f"9nt={pct_9nt:.2f}%, 5nt-total={pct_5nt_total:.2f}%, gap={diff:.2f}%. "
            f"Barcode present and correctly placed in {pct_5nt_total:.2f}% of reads; "
            f"TCGT suffix degraded/mutated in ~{diff:.1f}% of those reads."
        )

    return " | ".join(notes_parts), trimming_rec


# ============================================================
# MAIN
# ============================================================

def main():
    start_time = time.time()

    print("=" * 90)
    print("GSE53080 -- 28 WARNING FILES: 5nt vs 9nt DIAGNOSTIC v4")
    print(f"Timestamp:     {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Base path:     {BASE_PATH}")
    print(f"Input report:  {REPORT_PATH}")
    print(f"Output report: {OUTPUT_PATH}")
    print("=" * 90)
    print("\nMethod: single-pass sliding search.")
    print("  9nt checked first. If not found, 5nt counted with seq.count().")
    print("  Reads with 5nt appearing 2+ times flagged as false positive contamination.")
    print("  Trimming: always at 5nt barcode position, regardless of conclusion.")
    print("=" * 90)

    base = Path(BASE_PATH)
    if not base.exists():
        print(f"\nERROR: Base path does not exist: {base}")
        sys.exit(1)

    if not REPORT_PATH.exists():
        print(f"\nERROR: Previous barcode verification report not found:")
        print(f"       {REPORT_PATH}")
        print(f"\nPlease run verify_barcodes.py first, then re-run this script.")
        sys.exit(1)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    warning_files = load_warning_files_from_report(REPORT_PATH)
    n_warn = len(warning_files)
    print(f"\nLoaded {n_warn} WARNING file(s) from {REPORT_PATH.name}")
    if n_warn != 28:
        print(f"  WARNING: Expected 28 WARNING files; found {n_warn}. Proceeding.")

    output_rows = []
    total_files = n_warn

    print(f"\n{'#':>3} {'SRR_ID':<12} {'Folder':<25} {'TotReads':>10} "
          f"{'PrevPct':>8} {'9nt%':>8} {'5ntTot%':>8} {'Gap%':>8} "
          f"{'DblHit%':>8} {'9ntPos':>8} {'5ntPos':>8} {'Conclusion'}")
    print("-" * 110)

    for idx, info in enumerate(warning_files, 1):
        srr_id       = info["srr_id"]
        folder_rel   = info["folder"]
        adapter_id   = info["assigned_adapter_id"]
        barcode_5nt  = info["assigned_barcode"]
        previous_pct = info["previous_pct"]

        # Validate adapter ID
        if not adapter_id or adapter_id not in ADAPTERS:
            print(f"{idx:>3} {srr_id:<12} {folder_rel:<25} INVALID ADAPTER ID: {adapter_id}")
            output_rows.append({
                "SRR_ID": srr_id, "Folder": folder_rel,
                "Assigned_Barcode": barcode_5nt, "Assigned_Adapter_ID": adapter_id,
                "Full_Adapter": "N/A", "Total_Reads": 0,
                "Reads_9nt": 0, "Reads_5nt_Only": 0, "Reads_5nt_Total": 0,
                "Double_Hit_Reads": 0, "Double_Hit_Pct": 0.0,
                "Pct_9nt_Previous": previous_pct, "Pct_9nt": 0.0,
                "Pct_5nt_Total": 0.0, "Pct_5nt_Only": 0.0,
                "Difference_5nt_minus_9nt": 0.0, "Status": "WARNING",
                "Notes": f"Invalid adapter ID: {adapter_id}",
                "Barcode_Position_9nt_Median": "N/A",
                "Barcode_Position_5nt_Only_Median": "N/A",
                "RNA_Insert_Length_9nt_Median": "N/A",
                "RNA_Insert_Length_5nt_Only_Median": "N/A",
                "Read_Length_Median": "N/A",
                "Trimming_Recommendation": "Investigate -- invalid adapter assignment",
            })
            continue

        sig_9nt      = ID_TO_SIG[adapter_id]
        full_adapter = ID_TO_FULL[adapter_id]

        filepath = find_fastq(base, folder_rel, srr_id)
        if filepath is None:
            print(f"{idx:>3} {srr_id:<12} {folder_rel:<25} FILE NOT FOUND")
            output_rows.append({
                "SRR_ID": srr_id, "Folder": folder_rel,
                "Assigned_Barcode": barcode_5nt, "Assigned_Adapter_ID": adapter_id,
                "Full_Adapter": full_adapter, "Total_Reads": 0,
                "Reads_9nt": 0, "Reads_5nt_Only": 0, "Reads_5nt_Total": 0,
                "Double_Hit_Reads": 0, "Double_Hit_Pct": 0.0,
                "Pct_9nt_Previous": previous_pct, "Pct_9nt": 0.0,
                "Pct_5nt_Total": 0.0, "Pct_5nt_Only": 0.0,
                "Difference_5nt_minus_9nt": 0.0, "Status": "WARNING",
                "Notes": "FASTQ file not found",
                "Barcode_Position_9nt_Median": "N/A",
                "Barcode_Position_5nt_Only_Median": "N/A",
                "RNA_Insert_Length_9nt_Median": "N/A",
                "RNA_Insert_Length_5nt_Only_Median": "N/A",
                "Read_Length_Median": "N/A",
                "Trimming_Recommendation": "Investigate -- file missing",
            })
            continue

        print(f"{idx:>3} {srr_id:<12} scanning...", end="", flush=True)

        (total_reads, match_9nt, pos_9nt_list,
         match_5nt_only, pos_5nt_only_list,
         double_hit_reads, read_lengths) = scan_file(filepath, barcode_5nt, sig_9nt)

        match_5nt_total = match_9nt + match_5nt_only

        pct_9nt        = (match_9nt       / total_reads * 100.0) if total_reads > 0 else 0.0
        pct_5nt_total  = (match_5nt_total / total_reads * 100.0) if total_reads > 0 else 0.0
        pct_5nt_only   = (match_5nt_only  / total_reads * 100.0) if total_reads > 0 else 0.0
        double_hit_pct = (double_hit_reads / total_reads * 100.0) if total_reads > 0 else 0.0
        diff           = pct_5nt_total - pct_9nt

        checkpoint_diff = abs(pct_9nt - previous_pct)
        checkpoint_pass = checkpoint_diff <= 1.0

        pos_9nt_median      = compute_median(pos_9nt_list)
        pos_5nt_only_median = compute_median(pos_5nt_only_list)
        read_length_median  = compute_median(read_lengths)

        insert_9nt_median      = (pos_9nt_median - 1)      if pos_9nt_median      != "N/A" else "N/A"
        insert_5nt_only_median = (pos_5nt_only_median - 1) if pos_5nt_only_median != "N/A" else "N/A"

        notes, trimming_rec = build_conclusion(
            checkpoint_pass, checkpoint_diff,
            pct_9nt, previous_pct,
            pct_5nt_total, pct_5nt_only,
            double_hit_reads, total_reads,
            pos_9nt_median, pos_5nt_only_median
        )

        folder_short = (folder_rel
                        .replace("PLASMA_51nt/", "P/")
                        .replace("SERUM_51nt/",  "S/")
                        .replace("MYOCARDIUM_51nt/", "M51/")
                        .replace("MYOCARDIUM_36nt/", "M36/")
                        .replace("OTHER_TISSUES/",   "OT/"))

        conclusion_short = ("TCGT_ERR"   if "TCGT ERROR"        in notes
                       else "FALSE_POS"  if "FALSE POSITIVE"    in notes
                       else "REPRO_FAIL" if "REPRODUCIBILITY"   in notes
                       else "UNKNOWN")

        print(f"\r{idx:>3} {srr_id:<12} {folder_short:<25} {total_reads:>10,} "
              f"{previous_pct:>8.2f} {pct_9nt:>8.2f} {pct_5nt_total:>8.2f} {diff:>8.2f} "
              f"{double_hit_pct:>8.2f} {pos_9nt_median:>8} {pos_5nt_only_median:>8} "
              f"{conclusion_short}")

        output_rows.append({
            "SRR_ID":                           srr_id,
            "Folder":                           folder_rel,
            "Assigned_Barcode":                 barcode_5nt,
            "Assigned_Adapter_ID":              adapter_id,
            "Full_Adapter":                     full_adapter,
            "Total_Reads":                      total_reads,
            "Reads_9nt":                        match_9nt,
            "Reads_5nt_Only":                   match_5nt_only,
            "Reads_5nt_Total":                  match_5nt_total,
            "Double_Hit_Reads":                 double_hit_reads,
            "Double_Hit_Pct":                   round(double_hit_pct, 4),
            "Pct_9nt_Previous":                 round(previous_pct, 2),
            "Pct_9nt":                          round(pct_9nt, 2),
            "Pct_5nt_Total":                    round(pct_5nt_total, 2),
            "Pct_5nt_Only":                     round(pct_5nt_only, 2),
            "Difference_5nt_minus_9nt":         round(diff, 2),
            "Status":                           "WARNING",
            "Notes":                            notes,
            "Barcode_Position_9nt_Median":      pos_9nt_median,
            "Barcode_Position_5nt_Only_Median": pos_5nt_only_median,
            "RNA_Insert_Length_9nt_Median":     insert_9nt_median,
            "RNA_Insert_Length_5nt_Only_Median": insert_5nt_only_median,
            "Read_Length_Median":               read_length_median,
            "Trimming_Recommendation":          trimming_rec,
        })

    # Write CSV
    csv_fields = [
        "SRR_ID", "Folder", "Assigned_Barcode", "Assigned_Adapter_ID", "Full_Adapter",
        "Total_Reads", "Reads_9nt", "Reads_5nt_Only", "Reads_5nt_Total",
        "Double_Hit_Reads", "Double_Hit_Pct",
        "Pct_9nt_Previous", "Pct_9nt", "Pct_5nt_Total", "Pct_5nt_Only",
        "Difference_5nt_minus_9nt", "Status", "Notes",
        "Barcode_Position_9nt_Median", "Barcode_Position_5nt_Only_Median",
        "RNA_Insert_Length_9nt_Median", "RNA_Insert_Length_5nt_Only_Median",
        "Read_Length_Median", "Trimming_Recommendation"
    ]

    with open(OUTPUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        writer.writerows(output_rows)

    elapsed     = time.time() - start_time
    elapsed_str = str(timedelta(seconds=int(elapsed)))

    repro_fail = sum(1 for r in output_rows if "REPRODUCIBILITY" in r["Notes"])
    tcgt_error = sum(1 for r in output_rows if "TCGT ERROR"      in r["Notes"])
    false_pos  = sum(1 for r in output_rows if "FALSE POSITIVE"  in r["Notes"])

    print("\n" + "=" * 90)
    print("SUMMARY")
    print("=" * 90)
    print(f"Files scanned:                  {total_files}")
    print(f"Checkpoint passed (<=1% drift): {total_files - repro_fail}")
    print(f"Reproducibility FAIL:           {repro_fail}")
    print(f"TCGT ERROR confirmed:           {tcgt_error}")
    print(f"FALSE POSITIVE (partial):       {false_pos}")
    print(f"Elapsed time:                   {elapsed_str}")
    print(f"\nOutput saved: {OUTPUT_PATH.resolve()}")
    print("=" * 90)
    print("\nKEY:")
    print("  TCGT ERROR    -- barcode in ~100% reads, TCGT suffix degraded.")
    print("  FALSE POSITIVE -- 5nt barcode appears 2+ times in >10% of reads.")
    print("  Trimming: always cut at 5nt barcode position for both cases.")
    print("=" * 90)

    if repro_fail > 0:
        print("\nWARNING: Reproducibility check failed -- review before trimming!")
        sys.exit(1)
    else:
        print("\nAll checkpoints passed. Diagnostic complete. Ready for trimming.")


if __name__ == "__main__":
    main()
