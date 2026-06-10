#!/usr/bin/env python3
"""
GSE53080 FastQC + MultiQC Pipeline
==================================
Runs FastQC on all 180 analysis samples, then MultiQC on segregated batches.

Batch segregation (from SraRunTable metadata):
  - Batch A: HiSeq 2000 (150 samples) → split into serum (18) vs non-serum (132)
  - Batch B: GA IIx (30 samples) → all together
  - 5 excluded samples skipped (fetal cardiac, read-length mismatch)

Outputs:
  qc/fastqc/hiseq2000/serum/      ← 18 serum FastQC reports
  qc/fastqc/hiseq2000/non_serum/  ← 132 non-serum FastQC reports
  qc/fastqc/gaiix/all/            ← 30 GA IIx FastQC reports
  qc/multiqc/batch_hiseq/         ← MultiQC on 132 non-serum HiSeq
  qc/multiqc/batch_hiseq_serum/   ← MultiQC on 18 serum HiSeq
  qc/multiqc/batch_gaiix/         ← MultiQC on 30 GA IIx

Usage:
  python run_fastqc_multiqc.py --base /home/saravana/Downloads/GSE53080

  # Or skip FastQC and only run MultiQC on existing reports:
  python run_fastqc_multiqc.py --base /home/saravana/Downloads/GSE53080 --skip-fastqc
"""

import argparse
import os
import subprocess
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
EXCLUDED_SRR = {
    "SRR1044510", "SRR1044511", "SRR1044512",
    "SRR1044516", "SRR1044517",
}

# SRR → batch/subfolder mapping (generated from SraRunTable.csv, 180 samples)
BATCH_MAP = {
  "SRR1044415": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044416": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044417": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044418": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044419": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044420": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044421": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044422": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044423": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044424": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044425": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044426": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044427": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044428": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044429": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044430": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044431": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044432": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044433": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044434": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044435": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044436": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044437": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044438": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044439": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044440": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044441": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044442": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044443": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044444": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044445": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044446": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044447": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044448": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044449": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044450": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044451": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044452": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044453": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044454": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044455": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044456": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044457": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044458": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044459": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044460": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044461": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044462": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044463": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044464": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044465": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044466": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044467": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044468": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044469": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044470": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044471": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044472": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044473": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044474": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044475": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044476": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044477": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044478": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044479": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044480": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044481": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044482": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044483": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044484": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044485": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044486": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044487": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044488": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044489": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044490": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044491": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044492": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044493": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044494": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044495": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044496": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044497": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044498": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044499": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044500": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044501": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044502": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044503": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044504": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044505": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044506": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044507": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044508": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044509": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044513": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044514": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044515": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044518": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044519": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044520": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044521": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044522": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044523": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044524": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044525": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044526": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044527": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044528": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044529": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044530": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044531": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044532": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044533": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044534": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044535": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044536": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044537": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044538": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044539": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044540": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044541": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044542": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044543": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044544": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044545": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044546": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044547": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044548": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044549": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044550": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044551": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044552": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044553": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044554": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044555": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044556": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044557": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044558": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044559": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044560": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044561": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044562": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044563": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044564": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044565": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044566": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044567": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044568": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044569": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044570": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044571": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044572": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044573": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044574": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044575": {
    "batch": "hiseq2000",
    "sub": "serum"
  },
  "SRR1044576": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044577": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044578": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044579": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044580": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044581": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044582": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044583": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044584": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044585": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044586": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044587": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044588": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044589": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044590": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044591": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044592": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044593": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044594": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044595": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044596": {
    "batch": "hiseq2000",
    "sub": "non_serum"
  },
  "SRR1044597": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044598": {
    "batch": "gaiix",
    "sub": "all"
  },
  "SRR1044599": {
    "batch": "gaiix",
    "sub": "all"
  }
}

# ---------------------------------------------------------------------------
# UTILITIES
# ---------------------------------------------------------------------------

def find_fastq_files(base_dir):
    """Walk base_dir and return dict {SRR_ID: absolute_path}."""
    base = Path(base_dir)
    fastq_map = {}
    for path in base.rglob("*.fastq.gz"):
        fname = path.name
        if fname.startswith("SRR"):
            srr = fname.split("_")[0].split(".")[0]
            fastq_map[srr] = str(path.resolve())
    return fastq_map


def run_fastqc(args):
    """Worker function for parallel FastQC."""
    srr, fq_path, out_dir = args
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "fastqc",
        "--outdir", str(out_dir),
        fq_path,
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return (srr, "OK", None)
    except subprocess.CalledProcessError as e:
        return (srr, "FAIL", e.stderr)
    except FileNotFoundError:
        return (srr, "MISSING", "fastqc not found in PATH")


def run_multiqc(input_dir, out_dir, title):
    """Run MultiQC on a directory of FastQC outputs."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "multiqc",
        str(input_dir),
        "--outdir", str(out_dir),
        "--filename", f"multiqc_report_{title}",
        "--title", title,
        "--force",
    ]
    print(f"\n[MultiQC] {title}")
    print(f"   Input : {input_dir}")
    print(f"   Output: {out_dir}")
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run FastQC + MultiQC on GSE53080 with batch segregation"
    )
    parser.add_argument(
        "--base", required=True,
        help="Base directory containing GSE53080 fastq files"
    )
    parser.add_argument(
        "--qc-out", default="qc",
        help="Output directory for QC reports (default: ./qc)"
    )
    parser.add_argument(
        "--threads", type=int, default=12,
        help="Number of parallel FastQC jobs (default: 8 or CPU count)"
    )
    parser.add_argument(
        "--skip-fastqc", action="store_true",
        help="Skip FastQC (only run MultiQC on existing reports)"
    )
    parser.add_argument(
        "--skip-multiqc", action="store_true",
        help="Skip MultiQC (only run FastQC)"
    )
    args = parser.parse_args()

    base_dir = Path(args.base)
    qc_root = Path(args.qc_out).resolve()

    if not base_dir.exists():
        print(f"[ERROR] Base directory not found: {base_dir}")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 1. Discover fastq files
    # ------------------------------------------------------------------
    print(f"[1/4] Scanning {base_dir} for fastq.gz files...")
    fastq_map = find_fastq_files(base_dir)
    print(f"       Found {len(fastq_map)} fastq.gz files")

    # ------------------------------------------------------------------
    # 2. Match to 180-sample analysis set
    # ------------------------------------------------------------------
    analysis_srr = set(BATCH_MAP.keys())
    found_srr = set(fastq_map.keys())

    matched = analysis_srr & found_srr
    missing = analysis_srr - found_srr
    extra = found_srr - analysis_srr - EXCLUDED_SRR

    print(f"\n[2/4] Analysis set: {len(analysis_srr)} samples")
    print(f"       Matched   : {len(matched)}")
    print(f"       Missing   : {len(missing)}")
    if missing:
        for s in sorted(missing):
            print(f"         - {s}")
    print(f"       Extra     : {len(extra)} (not in 180-set, not excluded)")
    if extra:
        for s in sorted(extra)[:10]:
            print(f"         - {s}")

    if len(matched) != 180:
        print(f"[WARN] Expected 180 matched files, found {len(matched)}. Continuing with matched set.")

    # ------------------------------------------------------------------
    # 3. Run FastQC (parallel)
    # ------------------------------------------------------------------
    if not args.skip_fastqc:
        print(f"\n[3/4] Running FastQC on {len(matched)} files (threads={args.threads})...")

        tasks = []
        for srr in sorted(matched):
            info = BATCH_MAP[srr]
            out_dir = qc_root / "fastqc" / info["batch"] / info["sub"]
            tasks.append((srr, fastq_map[srr], out_dir))

        with Pool(args.threads) as pool:
            results = pool.map(run_fastqc, tasks)

        ok = sum(1 for _, status, _ in results if status == "OK")
        fail = sum(1 for _, status, _ in results if status == "FAIL")
        miss = sum(1 for _, status, _ in results if status == "MISSING")
        print(f"       Done: {ok} OK, {fail} FAILED, {miss} MISSING fastqc binary")

        if fail > 0:
            print("\nFailed samples:")
            for srr, status, err in results:
                if status == "FAIL":
                    print(f"  {srr}: {err}")
    else:
        print("\n[3/4] Skipping FastQC (--skip-fastqc)")

    # ------------------------------------------------------------------
    # 4. Run MultiQC on segregated batches
    # ------------------------------------------------------------------
    if not args.skip_multiqc:
        print("\n[4/4] Running MultiQC on segregated batches...")

        # Batch A1: HiSeq 2000 non-serum (132 samples)
        run_multiqc(
            qc_root / "fastqc" / "hiseq2000" / "non_serum",
            qc_root / "multiqc" / "batch_hiseq",
            "GSE53080_HiSeq2000_NonSerum",
        )

        # Batch A2: HiSeq 2000 serum (18 samples)
        run_multiqc(
            qc_root / "fastqc" / "hiseq2000" / "serum",
            qc_root / "multiqc" / "batch_hiseq_serum",
            "GSE53080_HiSeq2000_Serum",
        )

        # Batch B: GA IIx all (30 samples)
        run_multiqc(
            qc_root / "fastqc" / "gaiix" / "all",
            qc_root / "multiqc" / "batch_gaiix",
            "GSE53080_GAIIx_All",
        )

        print("\n[✓] All MultiQC reports generated.")
    else:
        print("\n[4/4] Skipping MultiQC (--skip-multiqc)")

    print(f"\nOutput structure:")
    print(f"  {qc_root}/fastqc/hiseq2000/serum/      ← 18 samples")
    print(f"  {qc_root}/fastqc/hiseq2000/non_serum/  ← 132 samples")
    print(f"  {qc_root}/fastqc/gaiix/all/            ← 30 samples")
    print(f"  {qc_root}/multiqc/batch_hiseq/         ← 132-sample report")
    print(f"  {qc_root}/multiqc/batch_hiseq_serum/   ← 18-sample report")
    print(f"  {qc_root}/multiqc/batch_gaiix/         ← 30-sample report")


if __name__ == "__main__":
    main()
