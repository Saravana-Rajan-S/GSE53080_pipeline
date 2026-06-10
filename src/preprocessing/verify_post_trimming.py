#!/usr/bin/env python3
"""
verify_post_trimming.py — Post-Trimming Verification (Step 6)
GSE53080 Pipeline — FINAL VERSION

Purpose:
  1. Scan all 720 trimmed files across 4 folders for residual 5nt barcodes
  2. Sliding window exact match (0 mismatch), same methodology as verify_barcodes.py
  3. Compute length distributions per file
  4. Compare 9nt match vs. 5nt rescue method side-by-side
  5. Identify rescued files and quantify recovery from 9nt-discarded reads

Author: Saravana Rajan S
Date: 2026-06-09
"""

import gzip
import csv
import statistics
import time
from pathlib import Path
from datetime import datetime, timedelta
from collections import Counter

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent.parent.parent  # repo root

# Original FASTQ files (35 subdirectories)
DOWNLOADS_PATH = Path("/home/saravana/Downloads/GSE53080")

# Pipeline paths
PIPELINE_PATH = Path("/home/saravana/Documents/GSE53080_pipeline")
DOCS_DIR = PIPELINE_PATH / "docs"
CONFIG_DIR = PIPELINE_PATH / "config"

# Trimmed output folders
TRIMMED_9NT_MATCH_DIR     = PIPELINE_PATH / "trimmed" / "9nt_match"
TRIMMED_RESCUED_INPUT_DIR = PIPELINE_PATH / "trimmed" / "rescued_input"
TRIMMED_5NT_RESCUED_DIR   = PIPELINE_PATH / "trimmed" / "5nt_rescued"
TRIMMED_5NT_MATCHED_DIR   = PIPELINE_PATH / "trimmed" / "5nt_matched"

# Input reference files
BARCODE_REPORT = DOCS_DIR / "barcode_verification_report.csv"
BARCODE_MAPPING = CONFIG_DIR / "barcode_mapping.csv"

# Output files
OUTPUT_VERIFICATION = DOCS_DIR / "post_trimming_verification.csv"
OUTPUT_COMPARISON = DOCS_DIR / "trimming_comparison_summary.csv"

# ──────────────────────────────────────────────────────────────────────────────
# 20 ADAPTERS: 5nt barcode + 9nt signature + full 26nt adapter
# ──────────────────────────────────────────────────────────────────────────────

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

# All 20 barcodes for cross-checking
ALL_BARCODES_5NT = [v["barcode"] for v in ADAPTERS.values()]

# ──────────────────────────────────────────────────────────────────────────────
# UTILITY FUNCTIONS
# ──────────────────────────────────────────────────────────────────────────────

def compute_median(values):
    """Compute median. Returns 0 if list is empty."""
    if not values:
        return 0
    return statistics.median(values)


def compute_mean(values):
    """Compute mean. Returns 0.0 if list is empty."""
    if not values:
        return 0.0
    return statistics.mean(values)


def pct(n, total):
    """Safe percentage: returns 0.0 if total is 0."""
    return round(n / total * 100, 4) if total > 0 else 0.0


def open_fastq(filepath):
    """Open FASTQ (gzip or plain text)."""
    if str(filepath).endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def scan_file_sliding(filepath, barcode_5nt):
    """
    SINGLE-PASS sliding window scan of a FASTQ file.

    Same methodology as verify_barcodes.py:
      - Sliding search: seq.find(barcode_5nt) anywhere in the read
      - Exact match, 0 mismatch
      - Counts reads containing the barcode at least once

    Returns dict with:
      total_reads, reads_with_barcode, barcode_pct,
      min_length, max_length, median_length, mean_length
    """
    total_reads = 0
    reads_with_barcode = 0
    lengths = []

    with open_fastq(filepath) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()

            # Guard against truncated files
            if not seq or not plus or not qual:
                break

            total_reads += 1
            lengths.append(len(seq))

            # Sliding search for exact 5nt barcode match
            if seq.find(barcode_5nt) != -1:
                reads_with_barcode += 1

    return {
        'total_reads': total_reads,
        'reads_with_barcode': reads_with_barcode,
        'barcode_pct': pct(reads_with_barcode, total_reads),
        'min_length': min(lengths) if lengths else 0,
        'max_length': max(lengths) if lengths else 0,
        'median_length': round(compute_median(lengths), 1),
        'mean_length': round(compute_mean(lengths), 2),
    }


def find_trimmed_file(directory, srr_id):
    """Find trimmed FASTQ file in a directory."""
    for suffix in ['.fastq.gz', '.fastq']:
        p = directory / f"{srr_id}{suffix}"
        if p.exists():
            return p
    return None


# ──────────────────────────────────────────────────────────────────────────────
# LOAD REFERENCE DATA
# ──────────────────────────────────────────────────────────────────────────────

def load_barcode_report(report_path):
    """
    Load barcode_verification_report.csv.
    Returns: {SRR_ID: {barcode, adapter_id, folder, total_reads, adapter_pct, status}}
    """
    data = {}
    with open(report_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            srr = row['SRR_ID'].strip()
            data[srr] = {
                'barcode': row['Assigned_Barcode'].strip(),
                'adapter_id': row['Assigned_Adapter_ID'].strip(),
                'folder': row['Folder'].strip(),
                'total_reads': int(row['Total_Reads']) if row['Total_Reads'] else 0,
                'adapter_pct': float(row['Adapter_Percentage']) if row['Adapter_Percentage'] else 0.0,
                'status': row['Status'].strip(),
            }
    return data


# ──────────────────────────────────────────────────────────────────────────────
# MAIN VERIFICATION
# ──────────────────────────────────────────────────────────────────────────────

def main():
    start_time = time.time()

    print("=" * 70)
    print("POST-TRIMMING VERIFICATION — GSE53080 Pipeline (Step 6)")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    print()
    print("Methodology: Sliding window exact match (0 mismatch) for 5nt barcode")
    print("Reference:   barcode_verification_report.csv")
    print("Folders:     9nt_match/, rescued_input/, 5nt_rescued/, 5nt_matched/")
    print("=" * 70)
    print()

    # Create output directories
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load reference data ──
    print("[1/4] Loading barcode verification report...")
    if not BARCODE_REPORT.exists():
        print(f"FATAL: Report not found: {BARCODE_REPORT}")
        print("Run verify_barcodes.py first.")
        return

    ref_data = load_barcode_report(BARCODE_REPORT)
    print(f"      Loaded {len(ref_data)} samples from reference report.")

    # ── Step 2: Verify trimmed directories ──
    print("\n[2/4] Checking trimmed output directories...")
    dirs_to_check = [
        ("9nt_match",     TRIMMED_9NT_MATCH_DIR),
        ("rescued_input", TRIMMED_RESCUED_INPUT_DIR),
        ("5nt_rescued",   TRIMMED_5NT_RESCUED_DIR),
        ("5nt_matched",   TRIMMED_5NT_MATCHED_DIR),
    ]
    for name, d in dirs_to_check:
        if d.exists():
            count = len(list(d.glob("SRR*.fastq*")))
            print(f"      {name:20s}: {d} ({count} files)")
        else:
            print(f"      {name:20s}: {d} — NOT FOUND")

    # ── Step 3: Process all samples ──
    print("\n[3/4] Scanning trimmed files...")
    verification_records = []
    comparison_records = []
    skipped = []

    total_srr = len(ref_data)
    sorted_srrs = sorted(ref_data.keys())

    for i, srr_id in enumerate(sorted_srrs, 1):
        if i % 20 == 0 or i == 1:
            print(f"      ... {i}/{total_srr} SRRs processed")

        ref = ref_data[srr_id]
        barcode_5nt = ref['barcode']
        adapter_id = ref['adapter_id']
        tissue_group = ref['folder']
        original_total = ref['total_reads']
        original_adapter_pct = ref['adapter_pct']
        original_status = ref['status']

        # Find files in each folder
        match_9nt_fp  = find_trimmed_file(TRIMMED_9NT_MATCH_DIR,     srr_id)
        rescued_fp    = find_trimmed_file(TRIMMED_RESCUED_INPUT_DIR,  srr_id)
        rescue_5nt_fp = find_trimmed_file(TRIMMED_5NT_RESCUED_DIR,    srr_id)
        matched_5nt_fp = find_trimmed_file(TRIMMED_5NT_MATCHED_DIR,   srr_id)

        # Check for missing critical files
        missing = []
        if not match_9nt_fp:  missing.append("9nt_match")
        if not rescued_fp:    missing.append("rescued_input")
        if not matched_5nt_fp: missing.append("5nt_matched")

        if missing:
            print(f"      SKIP {srr_id}: missing {', '.join(missing)}")
            skipped.append((srr_id, missing))
            continue

        # Scan each file (single-pass sliding search)
        match_9nt_stats = scan_file_sliding(match_9nt_fp, barcode_5nt)
        rescued_stats   = scan_file_sliding(rescued_fp, barcode_5nt)
        matched_5nt_stats = scan_file_sliding(matched_5nt_fp, barcode_5nt)

        # 5nt_rescued may not exist for all samples (only those needing rescue)
        rescue_5nt_stats = None
        if rescue_5nt_fp:
            rescue_5nt_stats = scan_file_sliding(rescue_5nt_fp, barcode_5nt)

        # ── Build verification records (one per file, consecutive by SRR) ──
        # Order: 9nt_match → rescued_input → 5nt_rescued (if exists) → 5nt_matched

        def make_verification_record(folder_name, file_path, stats, note):
            return {
                'SRR_ID': srr_id,
                'Folder': folder_name,
                'File_Name': file_path.name if file_path else 'N/A',
                'Assigned_Barcode_5nt': barcode_5nt,
                'Adapter_ID': adapter_id,
                'Tissue_Group': tissue_group,
                'Total_Reads': stats['total_reads'],
                'Reads_With_5nt_Barcode': stats['reads_with_barcode'],
                'Barcode_5nt_Pct': stats['barcode_pct'],
                'Min_Length': stats['min_length'],
                'Max_Length': stats['max_length'],
                'Median_Length': stats['median_length'],
                'Mean_Length': stats['mean_length'],
                'Original_Total_Reads': original_total,
                'Original_9nt_Adapter_Pct': original_adapter_pct,
                'Original_Status': original_status,
                'Notes': note,
            }

        # 1. 9nt_match/
        verification_records.append(make_verification_record(
            "9nt_match", match_9nt_fp, match_9nt_stats, "9nt exact match trimmed output"
        ))

        # 2. rescued_input/
        verification_records.append(make_verification_record(
            "rescued_input", rescued_fp, rescued_stats, "Leftover reads not matched by 9nt"
        ))

        # 3. 5nt_rescued/ (only if exists)
        if rescue_5nt_fp and rescue_5nt_stats:
            verification_records.append(make_verification_record(
                "5nt_rescued", rescue_5nt_fp, rescue_5nt_stats, "5nt rescue on leftovers"
            ))

        # 4. 5nt_matched/
        verification_records.append(make_verification_record(
            "5nt_matched", matched_5nt_fp, matched_5nt_stats, "Final merged output (9nt_match + 5nt_rescued)"
        ))

        # ── Build comparison record (one per SRR) ──
        rescued_5nt_reads = rescue_5nt_stats['total_reads'] if rescue_5nt_stats else 0
        delta_reads = matched_5nt_stats['total_reads'] - match_9nt_stats['total_reads']

        comparison_records.append({
            'SRR_ID': srr_id,
            'Assigned_Barcode_5nt': barcode_5nt,
            'Adapter_ID': adapter_id,
            'Tissue_Group': tissue_group,
            'Original_Status': original_status,
            'Original_Total_Reads': original_total,
            'Original_9nt_Adapter_Pct': original_adapter_pct,
            '9nt_Match_Reads': match_9nt_stats['total_reads'],
            '9nt_Match_Retention_Pct': pct(match_9nt_stats['total_reads'], original_total),
            'Rescued_Input_Reads': rescued_stats['total_reads'],
            'Rescued_Input_Pct': pct(rescued_stats['total_reads'], original_total),
            '5nt_Rescued_Reads': rescued_5nt_reads,
            '5nt_Rescued_Pct': pct(rescued_5nt_reads, original_total),
            '5nt_Matched_Reads': matched_5nt_stats['total_reads'],
            '5nt_Matched_Retention_Pct': pct(matched_5nt_stats['total_reads'], original_total),
            'Delta_Reads': delta_reads,
            'Delta_Pct': pct(delta_reads, original_total),
            'Rescue_Needed': 'YES' if rescued_5nt_reads > 0 else 'NO',
            '9nt_Match_Median_Length': match_9nt_stats['median_length'],
            '5nt_Matched_Median_Length': matched_5nt_stats['median_length'],
            '9nt_Match_Mean_Length': match_9nt_stats['mean_length'],
            '5nt_Matched_Mean_Length': matched_5nt_stats['mean_length'],
            'Residual_5nt_In_9nt_Match': match_9nt_stats['reads_with_barcode'],
            'Residual_5nt_In_5nt_Matched': matched_5nt_stats['reads_with_barcode'],
            'Residual_5nt_Pct_9nt_Match': match_9nt_stats['barcode_pct'],
            'Residual_5nt_Pct_5nt_Matched': matched_5nt_stats['barcode_pct'],
        })

    # ── Step 4: Write CSVs ──
    print("\n[4/4] Writing output files...")

    # Verification CSV (720 rows, consecutive by SRR)
    if verification_records:
        with open(OUTPUT_VERIFICATION, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(verification_records[0].keys()))
            writer.writeheader()
            writer.writerows(verification_records)
        print(f"      ✓ {OUTPUT_VERIFICATION}")
        print(f"        ({len(verification_records)} rows — consecutive SRR ordering, up to 4 rows per SRR)")

    # Comparison CSV (180 rows)
    if comparison_records:
        with open(OUTPUT_COMPARISON, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(comparison_records[0].keys()))
            writer.writeheader()
            writer.writerows(comparison_records)
        print(f"      ✓ {OUTPUT_COMPARISON}")
        print(f"        ({len(comparison_records)} rows — side-by-side comparison)")

    # ── Summary ──
    elapsed = time.time() - start_time
    elapsed_str = str(timedelta(seconds=int(elapsed)))

    print(f"\n{'=' * 70}")
    print("VERIFICATION SUMMARY")
    print(f"{'=' * 70}")
    print(f"SRRs processed     : {len(comparison_records)}/180")
    if skipped:
        print(f"Skipped            : {len(skipped)} SRRs")
        for srr, missing in skipped:
            print(f"    {srr}: missing {', '.join(missing)}")
    else:
        print(f"Skipped            : 0")

    # Rescue statistics
    rescued_srrs = [r for r in comparison_records if r['Rescue_Needed'] == 'YES']
    total_rescued = sum(r['5nt_Rescued_Reads'] for r in rescued_srrs)

    print(f"\n5nt RESCUE STATISTICS:")
    print(f"  Files with rescue  : {len(rescued_srrs)}")
    if rescued_srrs:
        print(f"  SRR IDs            : {[r['SRR_ID'] for r in rescued_srrs]}")
    print(f"  Total rescued reads: {total_rescued:,}")
    if rescued_srrs:
        avg_rescue_pct = statistics.mean([r['5nt_Rescued_Pct'] for r in rescued_srrs])
        print(f"  Avg rescue %       : {avg_rescue_pct:.2f}%")

    # Residual barcode check
    match_9nt_total = sum(r['9nt_Match_Reads'] for r in comparison_records)
    matched_5nt_total = sum(r['5nt_Matched_Reads'] for r in comparison_records)
    residual_9nt = sum(r['Residual_5nt_In_9nt_Match'] for r in comparison_records)
    residual_5nt = sum(r['Residual_5nt_In_5nt_Matched'] for r in comparison_records)

    print(f"\nRESIDUAL 5nt BARCODE CHECK (sliding search, exact match):")
    print(f"  9nt match strategy : {residual_9nt:,} reads ({residual_9nt/match_9nt_total*100:.4f}%)")
    print(f"  5nt matched output : {residual_5nt:,} reads ({residual_5nt/matched_5nt_total*100:.4f}%)")
    print(f"  (Should both be ~0 — any >0.1% indicates trimming bug)")

    # Length comparison
    match_9nt_medians   = [r['9nt_Match_Median_Length'] for r in comparison_records]
    matched_5nt_medians = [r['5nt_Matched_Median_Length'] for r in comparison_records]

    print(f"\nREAD LENGTH COMPARISON:")
    print(f"  9nt match median   : {statistics.median(match_9nt_medians):.1f} nt")
    print(f"  5nt matched median : {statistics.median(matched_5nt_medians):.1f} nt")
    print(f"  Delta              : {statistics.median(matched_5nt_medians) - statistics.median(match_9nt_medians):+.1f} nt")

    # 9nt missed reads summary
    total_delta = sum(r['Delta_Reads'] for r in comparison_records)
    total_original = sum(r['Original_Total_Reads'] for r in comparison_records)

    print(f"\n9nt MISSED READS (what 5nt rescue recovered):")
    print(f"  Total delta reads  : {total_delta:,}")
    print(f"  % of all original  : {total_delta/total_original*100:.2f}%")
    print(f"  (These are reads the 9nt pass discarded that 5nt rescue recovered)")

    print(f"\n{'=' * 70}")
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Elapsed: {elapsed_str}")
    print(f"{'=' * 70}")
    print("\nNext steps:")
    print("  1. Review post_trimming_verification.csv for per-file details")
    print("  2. Review trimming_comparison_summary.csv for side-by-side comparison")
    print("  3. Proceed to Step 7: Alignment")


if __name__ == "__main__":
    main()
