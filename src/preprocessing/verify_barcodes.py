#!/usr/bin/env python3
"""
GSE53080 Barcode Verification Script
=====================================
Location in repo: src/preprocessing/verify_barcodes.py

Scans ALL reads in each FASTQ file to identify the 5' adapter barcode.
The adapter is located at the 3' end of each read (after the RNA insert).
Performs a sliding search through the full read for the 9nt signature:
5nt sample-specific barcode + first 4nt of invariant Illumina adapter (TCGT).
Requires EXACT MATCH (0 mismatch) for confident assignment.

Read structure: [RNA insert ~18-30nt] + [5nt barcode] + TCGTATGCCGTCTTCTGCTTG

Outputs:
  1. docs/barcode_verification_report.csv  — Detailed per-file report with flags
  2. config/barcode_mapping.csv            — Clean SRR→Barcode mapping for trimming

Usage:
    cd ~/Documents/GSE53080_pipeline
    python3 src/preprocessing/verify_barcodes.py

    # With custom base path:
    BASE_PATH=/path/to/data python3 src/preprocessing/verify_barcodes.py
"""

import os
import sys
import gzip
import csv
import time
from pathlib import Path
from collections import Counter
from datetime import datetime, timedelta

# ============================================================
# CONFIGURATION — anchored to script location for reliability
# ============================================================
SCRIPT_DIR = Path(__file__).resolve().parent.parent.parent  # repo root
BASE_PATH = os.environ.get("BASE_PATH", "/home/saravana/Downloads/GSE53080")
REPORT_PATH = SCRIPT_DIR / "docs" / "barcode_verification_report.csv"
MAPPING_PATH = SCRIPT_DIR / "config" / "barcode_mapping.csv"

# 20 adapters: 5nt barcode + full 26nt adapter sequence
ADAPTERS = {
    "Ad01":  {"barcode": "TCACT", "signature": "TCACTTCGT", "full": "TCACTTCGTATGCCGTCTTCTGCTTG"},
    "Ad02":  {"barcode": "TCATC", "signature": "TCATCTCGT", "full": "TCATCTCGTATGCCGTCTTCTGCTTG"},
    "Ad03":  {"barcode": "TCCAC", "signature": "TCCACTCGT", "full": "TCCACTCGTATGCCGTCTTCTGCTTG"},
    "Ad04":  {"barcode": "TCCGT", "signature": "TCCGTTCGT", "full": "TCCGTTCGTATGCCGTCTTCTGCTTG"},
    "Ad05":  {"barcode": "TCCTA", "signature": "TCCTATCGT", "full": "TCCTATCGTATGCCGTCTTCTGCTTG"},
    "Ad06":  {"barcode": "TCGAT", "signature": "TCGATTCGT", "full": "TCGATTCGTATGCCGTCTTCTGCTTG"},
    "Ad07":  {"barcode": "TCGCG", "signature": "TCGCGTCGT", "full": "TCGCGTCGTATGCCGTCTTCTGCTTG"},
    "Ad08":  {"barcode": "TCTAG", "signature": "TCTAGTCGT", "full": "TCTAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad09":  {"barcode": "TCTCC", "signature": "TCTCCTCGT", "full": "TCTCCTCGTATGCCGTCTTCTGCTTG"},
    "Ad10":  {"barcode": "TCTGA", "signature": "TCTGATCGT", "full": "TCTGATCGTATGCCGTCTTCTGCTTG"},
    "Ad11":  {"barcode": "TTAAG", "signature": "TTAAGTCGT", "full": "TTAAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad12":  {"barcode": "TAACG", "signature": "TAACGTCGT", "full": "TAACGTCGTATGCCGTCTTCTGCTTG"},
    "Ad13":  {"barcode": "TAATA", "signature": "TAATATCGT", "full": "TAATATCGTATGCCGTCTTCTGCTTG"},
    "Ad14":  {"barcode": "TAGAG", "signature": "TAGAGTCGT", "full": "TAGAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad15":  {"barcode": "TAGGA", "signature": "TAGGATCGT", "full": "TAGGATCGTATGCCGTCTTCTGCTTG"},
    "Ad16":  {"barcode": "TATCA", "signature": "TATCATCGT", "full": "TATCATCGTATGCCGTCTTCTGCTTG"},
    "Ad17":  {"barcode": "TGATG", "signature": "TGATGTCGT", "full": "TGATGTCGTATGCCGTCTTCTGCTTG"},
    "Ad18":  {"barcode": "TGTGT", "signature": "TGTGTTCGT", "full": "TGTGTTCGTATGCCGTCTTCTGCTTG"},
    "Ad19":  {"barcode": "TTACA", "signature": "TTACATCGT", "full": "TTACATCGTATGCCGTCTTCTGCTTG"},
    "Ad20":  {"barcode": "TTGGT", "signature": "TTGGTTCGT", "full": "TTGGTTCGTATGCCGTCTTCTGCTTG"},
}

# Build reverse lookup: signature -> adapter ID
SIGNATURE_TO_ID = {v["signature"]: k for k, v in ADAPTERS.items()}
ALL_SIGNATURES = set(SIGNATURE_TO_ID.keys())

# ============================================================
# EXPECTED FOLDER STRUCTURE (from verify_gse53080.py)
# ============================================================
EXPECTED_STRUCTURE = {
    "MYOCARDIUM_51nt/01_NF": ["SRR1044415", "SRR1044416", "SRR1044417", "SRR1044418", "SRR1044419", "SRR1044420", "SRR1044421"],
    "MYOCARDIUM_51nt/02_FETAL": ["SRR1044513", "SRR1044515"],
    "MYOCARDIUM_51nt/03_DCM_ADVANCED": ["SRR1044452", "SRR1044454", "SRR1044457", "SRR1044458", "SRR1044461", "SRR1044469", "SRR1044474", "SRR1044475", "SRR1044476", "SRR1044477", "SRR1044478", "SRR1044479", "SRR1044482", "SRR1044485", "SRR1044487", "SRR1044489", "SRR1044494", "SRR1044495", "SRR1044496", "SRR1044500", "SRR1044502"],
    "MYOCARDIUM_51nt/04_DCM_EXPLANT": ["SRR1044451", "SRR1044453", "SRR1044460", "SRR1044468", "SRR1044488", "SRR1044499", "SRR1044501"],
    "MYOCARDIUM_51nt/05_ICM_ADVANCED": ["SRR1044527", "SRR1044524", "SRR1044529", "SRR1044535", "SRR1044536", "SRR1044544", "SRR1044549", "SRR1044551", "SRR1044555", "SRR1044560", "SRR1044566", "SRR1044577", "SRR1044578"],
    "MYOCARDIUM_51nt/06_ICM_EXPLANT": ["SRR1044523", "SRR1044528", "SRR1044534", "SRR1044550", "SRR1044565", "SRR1044576"],
    "MYOCARDIUM_36nt/01_NF": ["SRR1044422", "SRR1044423", "SRR1044424"],
    "MYOCARDIUM_36nt/02_DCM_ADVANCED": ["SRR1044456"],
    "MYOCARDIUM_36nt/03_DCM_EXPLANT": ["SRR1044455"],
    "MYOCARDIUM_36nt/04_ICM_ADVANCED": ["SRR1044522"],
    "MYOCARDIUM_36nt/05_ICM_EXPLANT": ["SRR1044521"],
    "PLASMA_51nt/01_NF": ["SRR1044425", "SRR1044427", "SRR1044429", "SRR1044431", "SRR1044433", "SRR1044434", "SRR1044435", "SRR1044436", "SRR1044437", "SRR1044438", "SRR1044439", "SRR1044440", "SRR1044441", "SRR1044442", "SRR1044443", "SRR1044444", "SRR1044445"],
    "PLASMA_51nt/02_DCM_ADVANCED": ["SRR1044464", "SRR1044465", "SRR1044471", "SRR1044481", "SRR1044484", "SRR1044486", "SRR1044491", "SRR1044498", "SRR1044505"],
    "PLASMA_51nt/03_DCM_3M_LVAD": ["SRR1044462", "SRR1044483", "SRR1044503"],
    "PLASMA_51nt/04_DCM_6M_LVAD": ["SRR1044480", "SRR1044497", "SRR1044504"],
    "PLASMA_51nt/05_DCM_EXPLANT": ["SRR1044463", "SRR1044470", "SRR1044490"],
    "PLASMA_51nt/06_ICM_ADVANCED": ["SRR1044526", "SRR1044531", "SRR1044539", "SRR1044540", "SRR1044546", "SRR1044548", "SRR1044554", "SRR1044557", "SRR1044559", "SRR1044562", "SRR1044564", "SRR1044568", "SRR1044569", "SRR1044573", "SRR1044580", "SRR1044582", "SRR1044584", "SRR1044586"],
    "PLASMA_51nt/07_ICM_3M_LVAD": ["SRR1044537", "SRR1044547", "SRR1044552", "SRR1044563", "SRR1044579", "SRR1044581", "SRR1044585"],
    "PLASMA_51nt/08_ICM_6M_LVAD": ["SRR1044525", "SRR1044545", "SRR1044553", "SRR1044556", "SRR1044558", "SRR1044561", "SRR1044583"],
    "PLASMA_51nt/09_ICM_EXPLANT": ["SRR1044530", "SRR1044538", "SRR1044567", "SRR1044572"],
    "PLASMA_51nt/10_STABLE_HF": ["SRR1044506", "SRR1044507", "SRR1044508", "SRR1044509", "SRR1044587", "SRR1044588", "SRR1044589", "SRR1044590", "SRR1044591", "SRR1044592", "SRR1044593", "SRR1044594", "SRR1044595", "SRR1044596"],
    "SERUM_51nt/01_NF": ["SRR1044426", "SRR1044428", "SRR1044430", "SRR1044432"],
    "SERUM_51nt/02_DCM_ADVANCED": ["SRR1044467", "SRR1044473", "SRR1044493"],
    "SERUM_51nt/03_DCM_EXPLANT": ["SRR1044466", "SRR1044472", "SRR1044492"],
    "SERUM_51nt/04_ICM_ADVANCED": ["SRR1044533", "SRR1044543", "SRR1044571", "SRR1044575"],
    "SERUM_51nt/05_ICM_EXPLANT": ["SRR1044532", "SRR1044542", "SRR1044570", "SRR1044574"],
    "OTHER_TISSUES/ADIPOSE": ["SRR1044446"],
    "OTHER_TISSUES/BRAIN": ["SRR1044447"],
    "OTHER_TISSUES/LIVER": ["SRR1044448"],
    "OTHER_TISSUES/SKIN": ["SRR1044449"],
    "OTHER_TISSUES/SPLEEN": ["SRR1044450"],
    "OTHER_TISSUES/SKELETAL_MUSCLE": ["SRR1044514", "SRR1044459", "SRR1044541", "SRR1044518"],
    "OTHER_TISSUES/HUVEC": ["SRR1044519", "SRR1044520"],
    "OTHER_TISSUES/PBMC": ["SRR1044597"],
    "OTHER_TISSUES/RBC": ["SRR1044598", "SRR1044599"],
}

# ============================================================
# FUNCTIONS
# ============================================================

def find_fastq(base_path, folder_rel, srr_id):
    """Find the FASTQ file for a given SRR in a folder."""
    folder = Path(base_path) / folder_rel
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        candidate = folder / f"{srr_id}{ext}"
        if candidate.exists():
            return candidate
    return None


def scan_file(filepath):
    """
    Scan entire FASTQ file. Perform sliding search for exact 9nt signature
    matches anywhere in the read. Collect barcode start positions (1-based)
    and read lengths.
    Returns: total_reads, barcode_counts (Counter), reads_with_adapter,
             barcode_positions (list, 1-based), read_lengths (list)
    """
    total_reads = 0
    barcode_counts = Counter()
    reads_with_adapter = 0
    barcode_positions = []  # 1-based start positions of matched signatures
    read_lengths = []       # full read lengths

    opener = gzip.open if str(filepath).endswith(".gz") else open

    with opener(filepath, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()

            # Guard against truncated/corrupted FASTQ files
            if not seq or not plus or not qual:
                print(f"    ⚠️  Truncated FASTQ detected at read {total_reads + 1}. Stopping scan.")
                break

            total_reads += 1
            read_lengths.append(len(seq))

            # Sliding search — check all 20 signatures anywhere in the read
            found = False
            for sig, adapter_id in SIGNATURE_TO_ID.items():
                pos = seq.find(sig)
                if pos != -1:
                    reads_with_adapter += 1
                    barcode_counts[adapter_id] += 1
                    barcode_positions.append(pos + 1)  # 1-based position
                    found = True
                    break  # Only one signature should match per read

    return total_reads, barcode_counts, reads_with_adapter, barcode_positions, read_lengths


def compute_median(values):
    """Compute median of a list of numeric values."""
    if not values:
        return "N/A"
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    mid = n // 2
    if n % 2 == 1:
        return sorted_vals[mid]
    else:
        return (sorted_vals[mid - 1] + sorted_vals[mid]) / 2


def assign_barcode(barcode_counts, total_reads, reads_with_adapter, barcode_positions, read_lengths):
    """
    Determine dominant barcode and assign status flags.
    Computes barcode position median (1-based), RNA insert length median,
    and read length median.
    Returns: (assigned_id, barcode_seq, full_adapter, status, notes,
              barcode_pos_median, insert_length_median, read_length_median)
    """
    pos_median = compute_median(barcode_positions)
    read_length_median = compute_median(read_lengths)

    if pos_median != "N/A":
        insert_length_median = pos_median - 1
    else:
        insert_length_median = "N/A"

    if total_reads == 0:
        return ("N/A", "N/A", "N/A", "CRITICAL",
                "Empty file or no valid reads found",
                pos_median, insert_length_median, read_length_median)

    adapter_pct = (reads_with_adapter / total_reads) * 100.0

    if reads_with_adapter == 0:
        return ("N/A", "N/A", "N/A", "CRITICAL",
                "No adapter detected (0%) — check file integrity",
                pos_median, insert_length_median, read_length_median)

    if not barcode_counts:
        return ("N/A", "N/A", "N/A", "CRITICAL",
                "No known barcode signatures found",
                pos_median, insert_length_median, read_length_median)

    # Find dominant barcode
    top_id, top_count = barcode_counts.most_common(1)[0]
    top_pct = (top_count / total_reads) * 100.0
    second_count = barcode_counts.most_common(2)[1][1] if len(barcode_counts) > 1 else 0

    barcode_seq = ADAPTERS[top_id]["barcode"]
    full_adapter = ADAPTERS[top_id]["full"]

    # Build notes list
    notes = []

    # Flag: low adapter rate
    if adapter_pct < 70.0:
        notes.append(f"Low adapter rate ({adapter_pct:.1f}%)")

    # Flag: near-universal adapter (possible adapter-dimer library)
    if adapter_pct > 99.9:
        notes.append("Near-universal adapter (>99.9%) — possible adapter-dimer library")

    # Flag: mixed barcodes (second barcode > 10% of top)
    if second_count > 0 and (second_count / top_count) > 0.10:
        second_id = barcode_counts.most_common(2)[1][0]
        notes.append(f"Mixed barcodes: {second_id}={second_count} vs {top_id}={top_count}")

    # Flag: weak dominant barcode
    if top_pct < 80.0:
        notes.append(f"Weak dominant barcode ({top_pct:.1f}%)")

    if not notes:
        status = "Clean"
        notes_str = "All checks passed"
    else:
        has_warning = any(
            "Low adapter rate" in n or "Mixed barcodes" in n or "Weak dominant" in n
            for n in notes
        )
        if has_warning:
            status = "WARNING"
        else:
            status = "INFO"
        notes_str = "; ".join(notes)

    return (top_id, barcode_seq, full_adapter, status, notes_str,
            pos_median, insert_length_median, read_length_median)


# ============================================================
# MAIN
# ============================================================

def main():
    start_time = time.time()

    print("=" * 70)
    print("GSE53080 BARCODE VERIFICATION")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Base path: {BASE_PATH}")
    print(f"Signature check: 9nt (5nt barcode + TCGT), sliding search through full read, EXACT MATCH (0 mismatch)")
    print(f"Output report:   {REPORT_PATH}")
    print(f"Output mapping:  {MAPPING_PATH}")
    print("=" * 70)

    base = Path(BASE_PATH)
    if not base.exists():
        print(f"\nERROR: Base path does not exist: {base}")
        print("Set BASE_PATH environment variable or edit the script.")
        sys.exit(1)

    # Prepare output directories
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    MAPPING_PATH.parent.mkdir(parents=True, exist_ok=True)

    report_rows = []
    mapping_rows = []

    total_files = sum(len(v) for v in EXPECTED_STRUCTURE.values())
    processed = 0

    for folder_rel, srr_list in EXPECTED_STRUCTURE.items():
        for srr_id in srr_list:
            processed += 1
            filepath = find_fastq(base, folder_rel, srr_id)

            if filepath is None:
                print(f"\n[{processed}/{total_files}] ❌ {srr_id} — FILE NOT FOUND in {folder_rel}")
                report_rows.append({
                    "SRR_ID": srr_id,
                    "Folder": folder_rel,
                    "Total_Reads": 0,
                    "Reads_With_Adapter": 0,
                    "Adapter_Percentage": 0.0,
                    "Assigned_Barcode": "N/A",
                    "Assigned_Adapter_ID": "N/A",
                    "Full_Adapter": "N/A",
                    "Status": "CRITICAL",
                    "Notes": "FASTQ file not found",
                    "Barcode_Start_Position_Median": "N/A",
                    "RNA_Insert_Length_Median": "N/A",
                    "Read_Length_Median": "N/A",
                })
                continue

            # Scan entire file
            print(f"[{processed}/{total_files}] Scanning {srr_id} ... ", end="", flush=True)
            total_reads, barcode_counts, reads_with_adapter, barcode_positions, read_lengths = scan_file(filepath)
            adapter_pct = (reads_with_adapter / total_reads * 100.0) if total_reads > 0 else 0.0

            # Assign barcode and flags
            (assigned_id, barcode_seq, full_adapter, status, notes,
             pos_median, insert_length_median, read_length_median) = assign_barcode(
                barcode_counts, total_reads, reads_with_adapter, barcode_positions, read_lengths
            )

            # Build position display string
            if pos_median != "N/A":
                pos_str = (f"PosMed={pos_median} InsMed={insert_length_median} "
                           f"ReadLenMed={read_length_median}")
            else:
                pos_str = "Pos=N/A"

            print(f"Total={total_reads:,}  Adapter={reads_with_adapter:,} ({adapter_pct:.1f}%)  "
                  f"Barcode={assigned_id}  {pos_str}  Status={status}")

            report_rows.append({
                "SRR_ID": srr_id,
                "Folder": folder_rel,
                "Total_Reads": total_reads,
                "Reads_With_Adapter": reads_with_adapter,
                "Adapter_Percentage": round(adapter_pct, 2),
                "Assigned_Barcode": barcode_seq,
                "Assigned_Adapter_ID": assigned_id,
                "Full_Adapter": full_adapter,
                "Status": status,
                "Notes": notes,
                "Barcode_Start_Position_Median": pos_median,
                "RNA_Insert_Length_Median": insert_length_median,
                "Read_Length_Median": read_length_median,
            })

            if status != "CRITICAL":
                mapping_rows.append({
                    "SRR_ID": srr_id,
                    "Barcode_9nt": ADAPTERS[assigned_id]["signature"],
                    "Adapter_ID": assigned_id,
                })

    # Write detailed report CSV
    report_fields = [
        "SRR_ID", "Folder", "Total_Reads", "Reads_With_Adapter",
        "Adapter_Percentage", "Assigned_Barcode", "Assigned_Adapter_ID",
        "Full_Adapter", "Status", "Notes",
        "Barcode_Start_Position_Median", "RNA_Insert_Length_Median", "Read_Length_Median"
    ]
    with open(REPORT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=report_fields)
        writer.writeheader()
        writer.writerows(report_rows)

    # Write clean mapping CSV for trimming
    mapping_fields = ["SRR_ID", "Barcode_9nt", "Adapter_ID"]
    with open(MAPPING_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=mapping_fields)
        writer.writeheader()
        writer.writerows(mapping_rows)

    # Summary
    clean_count = sum(1 for r in report_rows if r["Status"] == "Clean")
    warning_count = sum(1 for r in report_rows if r["Status"] == "WARNING")
    critical_count = sum(1 for r in report_rows if r["Status"] == "CRITICAL")
    info_count = sum(1 for r in report_rows if r["Status"] == "INFO")

    elapsed = time.time() - start_time
    elapsed_str = str(timedelta(seconds=int(elapsed)))

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total files scanned:      {total_files}")
    print(f"Clean:                    {clean_count}")
    print(f"INFO (minor note):        {info_count}")
    print(f"WARNING:                  {warning_count}")
    print(f"CRITICAL:                 {critical_count}")
    print(f"Elapsed time:             {elapsed_str}")
    print(f"\nReports saved:")
    print(f"  Detailed report:        {REPORT_PATH.resolve()}")
    print(f"  Barcode mapping:        {MAPPING_PATH.resolve()}")
    print("=" * 70)

    if critical_count > 0:
        print("\n⚠️  CRITICAL issues found — review report before trimming!")
        sys.exit(1)
    else:
        print("\n✅ All files scanned successfully. Ready for trimming.")


if __name__ == "__main__":
    main()
