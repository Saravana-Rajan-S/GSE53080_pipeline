#!/usr/bin/env python3
"""
GSE53080 Dual Trimming Pipeline — Sequential Two-Pass Workflow
===============================================================

Pass 1 (9nt match): 9nt exact match → trimmed/9nt_match/ + rescued_input/
Pass 2 (5nt match): 5nt exact match on leftovers → trimmed/5nt_rescued/
Pass 3 (Merge):     cat 9nt_match/ + 5nt_rescued/ → trimmed/5nt_matched/

Strategy A (9nt style):
  cutadapt -a [barcode]TCGT -O 9 -e 0 -m 15 --discard-untrimmed
  → --untrimmed-output saves TCGT-degraded reads for rescue

Strategy B (5nt rescue):
  cutadapt -a [barcode] -O 5 -e 0 -m 15 --discard-untrimmed
  → runs on rescued_input/ leftovers only

Input:
  - config/barcode_mapping.csv (SRR → barcode → adapter)
  - /home/saravana/Downloads/GSE53080/ (FASTQ files, recursive search)

Output:
  - trimmed/9nt_match/      ← 9nt exact only (180 files)
  - trimmed/rescued_input/  ← untrimmed leftovers from 9nt (180 files)
  - trimmed/5nt_rescued/    ← 5nt rescue on leftovers (180 files)
  - trimmed/5nt_matched/    ← final combined output (180 files)

Usage:
  # Full run (all 180 samples)
  python3 src/preprocessing/run_dual_trimming.py \
      --base /home/saravana/Downloads/GSE53080 \
      --mapping config/barcode_mapping.csv \
      --output trimmed/ \
      --threads 4

  # Dry run (print commands, don't execute)
  python3 src/preprocessing/run_dual_trimming.py \
      --base /home/saravana/Downloads/GSE53080 \
      --mapping config/barcode_mapping.csv \
      --output trimmed/ \
      --dry-run

  # Test on first 3 samples only
  python3 src/preprocessing/run_dual_trimming.py \
      --base /home/saravana/Downloads/GSE53080 \
      --mapping config/barcode_mapping.csv \
      --output trimmed/ \
      --max-samples 3
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from collections import defaultdict

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
ADAPTERS = {
    "Ad01": {"barcode": "TCACT", "signature_9nt": "TCACTTCGT", "full": "TCACTTCGTATGCCGTCTTCTGCTTG"},
    "Ad02": {"barcode": "TCATC", "signature_9nt": "TCATCTCGT", "full": "TCATCTCGTATGCCGTCTTCTGCTTG"},
    "Ad03": {"barcode": "TCCAC", "signature_9nt": "TCCACTCGT", "full": "TCCACTCGTATGCCGTCTTCTGCTTG"},
    "Ad04": {"barcode": "TCCGT", "signature_9nt": "TCCGTTCGT", "full": "TCCGTTCGTATGCCGTCTTCTGCTTG"},
    "Ad05": {"barcode": "TCCTA", "signature_9nt": "TCCTATCGT", "full": "TCCTATCGTATGCCGTCTTCTGCTTG"},
    "Ad06": {"barcode": "TCGAT", "signature_9nt": "TCGATTCGT", "full": "TCGATTCGTATGCCGTCTTCTGCTTG"},
    "Ad07": {"barcode": "TCGCG", "signature_9nt": "TCGCGTCGT", "full": "TCGCGTCGTATGCCGTCTTCTGCTTG"},
    "Ad08": {"barcode": "TCTAG", "signature_9nt": "TCTAGTCGT", "full": "TCTAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad09": {"barcode": "TCTCC", "signature_9nt": "TCTCCTCGT", "full": "TCTCCTCGTATGCCGTCTTCTGCTTG"},
    "Ad10": {"barcode": "TCTGA", "signature_9nt": "TCTGATCGT", "full": "TCTGATCGTATGCCGTCTTCTGCTTG"},
    "Ad11": {"barcode": "TTAAG", "signature_9nt": "TTAAGTCGT", "full": "TTAAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad12": {"barcode": "TAACG", "signature_9nt": "TAACGTCGT", "full": "TAACGTCGTATGCCGTCTTCTGCTTG"},
    "Ad13": {"barcode": "TAATA", "signature_9nt": "TAATATCGT", "full": "TAATATCGTATGCCGTCTTCTGCTTG"},
    "Ad14": {"barcode": "TAGAG", "signature_9nt": "TAGAGTCGT", "full": "TAGAGTCGTATGCCGTCTTCTGCTTG"},
    "Ad15": {"barcode": "TAGGA", "signature_9nt": "TAGGATCGT", "full": "TAGGATCGTATGCCGTCTTCTGCTTG"},
    "Ad16": {"barcode": "TATCA", "signature_9nt": "TATCATCGT", "full": "TATCATCGTATGCCGTCTTCTGCTTG"},
    "Ad17": {"barcode": "TGATG", "signature_9nt": "TGATGTCGT", "full": "TGATGTCGTATGCCGTCTTCTGCTTG"},
    "Ad18": {"barcode": "TGTGT", "signature_9nt": "TGTGTTCGT", "full": "TGTGTTCGTATGCCGTCTTCTGCTTG"},
    "Ad19": {"barcode": "TTACA", "signature_9nt": "TTACATCGT", "full": "TTACATCGTATGCCGTCTTCTGCTTG"},
    "Ad20": {"barcode": "TTGGT", "signature_9nt": "TTGGTTCGT", "full": "TTGGTTCGTATGCCGTCTTCTGCTTG"},
}

# ---------------------------------------------------------------------------
# UTILITIES
# ---------------------------------------------------------------------------
def find_fastq_file(base_dir, srr_id):
    """Find FASTQ file for SRR_ID recursively under base_dir."""
    base = Path(base_dir)
    for pattern in [f"{srr_id}*.fastq.gz", f"{srr_id}*.fq.gz"]:
        matches = list(base.rglob(pattern))
        if matches:
            return str(matches[0])
    return None


def run_cutadapt(cmd, dry_run=False):
    """Execute or print a cutadapt command."""
    if dry_run:
        print(f"  [DRY-RUN] {' '.join(cmd)}")
        return (0, "", "")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return (0, result.stdout, result.stderr)
    except subprocess.CalledProcessError as e:
        return (e.returncode, e.stdout, e.stderr)


# ---------------------------------------------------------------------------
# SINGLE-SAMPLE TRIMMING
# ---------------------------------------------------------------------------
def trim_sample(args):
    """
    Worker function for parallel trimming.
    args = (srr_id, adapter_id, fq_path, output_base, dry_run)
    """
    srr_id, adapter_id, fq_path, output_base, dry_run = args

    adapter = ADAPTERS[adapter_id]
    barcode_5nt = adapter["barcode"]
    signature_9nt = adapter["signature_9nt"]

    out_base = Path(output_base)
    match_9nt_dir = out_base / "9nt_match"
    rescued_dir = out_base / "rescued_input"
    rescue_5nt_dir = out_base / "5nt_rescued"
    matched_5nt_dir = out_base / "5nt_matched"

    # (directories created in main() before Pool)

    # File paths
    match_9nt_out = match_9nt_dir / f"{srr_id}.fastq.gz"
    rescued_out = rescued_dir / f"{srr_id}.fastq.gz"
    rescue_5nt_out = rescue_5nt_dir / f"{srr_id}.fastq.gz"
    matched_5nt_out = matched_5nt_dir / f"{srr_id}.fastq.gz"

    log_lines = []

    # =====================================================================
    # PASS 1: 9nt exact match
    # =====================================================================
    log_lines.append(f"[{srr_id}] Pass 1: 9nt exact match")

    cmd1 = [
        "cutadapt",
        "-a", signature_9nt,
        "-O", "9",
        "-e", "0",
        "-m", "15",
        "--untrimmed-output", str(rescued_out),
        "-o", str(match_9nt_out),
        fq_path,
    ]

    ret1, stdout1, stderr1 = run_cutadapt(cmd1, dry_run)
    if ret1 != 0:
        log_lines.append(f"  ERROR: Pass 1 failed: {stderr1}")
        return (srr_id, "FAIL_PASS1", "\n".join(log_lines))

    # Count reads in 9nt_match output (for reporting)
    if not dry_run:
        match_9nt_count = count_reads_fastq(match_9nt_out)
        rescued_count = count_reads_fastq(rescued_out)
        log_lines.append(f"  9nt matched:    {match_9nt_count} reads")
        log_lines.append(f"  Rescued input:  {rescued_count} reads")

    # =====================================================================
    # PASS 2: 5nt rescue (5nt exact on leftovers)
    # =====================================================================
    log_lines.append(f"[{srr_id}] Pass 2: 5nt rescue")

    # Only run if rescued_input has reads
    if not dry_run and rescued_count == 0:
        log_lines.append(f"  No rescued reads — skipping Pass 2")
        # Just copy 9nt_match output to 5nt_matched
        shutil.copy2(str(match_9nt_out), str(matched_5nt_out))
        return (srr_id, "OK_NO_RESCUE", "\n".join(log_lines))

    cmd2 = [
        "cutadapt",
        "-a", barcode_5nt,
        "-O", "5",
        "-e", "0",
        "-m", "15",
        "--discard-untrimmed",
        "-o", str(rescue_5nt_out),
        str(rescued_out),
    ]

    ret2, stdout2, stderr2 = run_cutadapt(cmd2, dry_run)
    if ret2 != 0:
        log_lines.append(f"  ERROR: Pass 2 failed: {stderr2}")
        return (srr_id, "FAIL_PASS2", "\n".join(log_lines))

    if not dry_run:
        rescue_5nt_count = count_reads_fastq(rescue_5nt_out)
        log_lines.append(f"  5nt rescued:    {rescue_5nt_count} reads")

    # =====================================================================
    # PASS 3: Merge 9nt_match + 5nt_rescued
    # =====================================================================
    log_lines.append(f"[{srr_id}] Pass 3: Merge")

    if dry_run:
        log_lines.append(f"  [DRY-RUN] cat {match_9nt_out} + {rescue_5nt_out} -> {matched_5nt_out}")
    else:
        # Merge using chunked copy (low memory)
        with open(matched_5nt_out, 'wb') as merged_fh:
            for src in [match_9nt_out, rescue_5nt_out]:
                if src.exists():
                    with open(src, 'rb') as src_fh:
                        shutil.copyfileobj(src_fh, merged_fh)

        matched_count = count_reads_fastq(matched_5nt_out)
        total_matched = match_9nt_count + (rescue_5nt_count if rescue_5nt_out.exists() else 0)
        log_lines.append(f"  Merged output:  {matched_count} reads")
        log_lines.append(f"  Total matched:  {total_matched} reads (9nt={match_9nt_count}, 5nt={rescue_5nt_count})")

    return (srr_id, "OK", "\n".join(log_lines))


def count_reads_fastq(filepath):
    """Count reads in a FASTQ file (handles .gz)."""
    import gzip
    opener = gzip.open if str(filepath).endswith('.gz') else open
    line_count = 0
    with opener(filepath, 'rt') as fh:
        for _ in fh:
            line_count += 1
    return line_count // 4


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="GSE53080 Dual Trimming — Sequential Two-Pass Workflow"
    )
    parser.add_argument(
        "--base", required=True,
        help="Base directory containing GSE53080 fastq files"
    )
    parser.add_argument(
        "--mapping", required=True,
        help="Path to barcode_mapping.csv"
    )
    parser.add_argument(
        "--output", default="trimmed",
        help="Output directory for trimmed files (default: trimmed/)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of parallel trimming jobs (default: 4)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing"
    )
    parser.add_argument(
        "--max-samples", type=int, default=None,
        help="Limit to first N samples for testing"
    )
    args = parser.parse_args()

    base_dir = Path(args.base)
    mapping_path = Path(args.mapping)
    output_base = Path(args.output)

    if not mapping_path.exists():
        print(f"[ERROR] Mapping file not found: {mapping_path}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 1. Load barcode mapping
    # ------------------------------------------------------------------
    samples = []
    with open(mapping_path, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            srr = row['SRR_ID']
            adapter_id = row.get('Adapter_ID', row.get('adapter_id', ''))

            # Find FASTQ file
            fq_path = find_fastq_file(base_dir, srr)
            if not fq_path:
                print(f"[WARN] No FASTQ found for {srr}")
                continue

            samples.append((srr, adapter_id, fq_path))

    if args.max_samples:
        samples = samples[:args.max_samples]

    print(f"[INFO] Loaded {len(samples)} samples from mapping")

    # ------------------------------------------------------------------
    # 2. Run trimming (parallel or sequential)
    # ------------------------------------------------------------------
    # Create output directories before spawning workers
    for subdir in ["9nt_match", "rescued_input", "5nt_rescued", "5nt_matched"]:
        (output_base / subdir).mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Starting trimming with {args.threads} threads...")
    print(f"[INFO] Output: {output_base.resolve()}")
    print()

    tasks = [(srr, aid, fq, str(output_base), args.dry_run) for srr, aid, fq in samples]

    if args.dry_run:
        # Sequential for dry-run (cleaner output)
        results = [trim_sample(t) for t in tasks]
    else:
        with Pool(args.threads) as pool:
            results = pool.map(trim_sample, tasks)

    # ------------------------------------------------------------------
    # 3. Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("TRIMMING SUMMARY")
    print("=" * 70)

    ok_count = sum(1 for _, status, _ in results if status == "OK")
    ok_no_rescue = sum(1 for _, status, _ in results if status == "OK_NO_RESCUE")
    fail_count = sum(1 for _, status, _ in results if status.startswith("FAIL"))

    print(f"  Total samples:    {len(results)}")
    print(f"  OK (with rescue): {ok_count}")
    print(f"  OK (no rescue):   {ok_no_rescue}")
    print(f"  Failed:           {fail_count}")
    print()

    if fail_count > 0:
        print("  Failed samples:")
        for srr, status, log in results:
            if status.startswith("FAIL"):
                print(f"    {srr}: {status}")
                print(f"      {log}")

    print(f"\nOutput structure:")
    print(f"  {output_base}/9nt_match/       ← 9nt exact only")
    print(f"  {output_base}/rescued_input/   ← untrimmed leftovers")
    print(f"  {output_base}/5nt_rescued/     ← 5nt rescue on leftovers")
    print(f"  {output_base}/5nt_matched/     ← final combined output")
    print("=" * 70)

    if not args.dry_run and fail_count == 0:
        print("\n✅ Trimming complete. Next: run post-trimming verification.")


if __name__ == "__main__":
    main()
