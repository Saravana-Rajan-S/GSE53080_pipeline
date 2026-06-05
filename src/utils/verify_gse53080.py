#!/usr/bin/env python3
"""
GSE53080 Folder Structure Verification Script
===============================================
Location in repo: src/utils/verify_gse53080.py

This script validates your folder structure against the metadata-verified
assignments from docs/folder_assignments.md.

Usage:
    python3 src/utils/verify_gse53080.py

    # Or from repo root:
    python3 verify_gse53080.py

    # With custom base path:
    BASE_PATH=/path/to/data python3 -c "import os; os.environ['BASE_PATH']='/path'; exec(open('verify_gse53080.py').read())"

"""

import os
import sys
from pathlib import Path
from datetime import datetime

# ============================================================
# CONFIGURATION — CHANGE THIS TO YOUR ACTUAL PATH
# ============================================================
BASE_PATH = "/home/saravana/Downloads/GSE53080"
# Alternative: BASE_PATH = os.path.expanduser("~/Downloads/GSE53080")

# ============================================================
# EXPECTED FOLDER STRUCTURE (metadata-verified)
# ============================================================
EXPECTED_STRUCTURE = {
    # MYOCARDIUM_51nt
    "MYOCARDIUM_51nt/01_NF": [
        "SRR1044415", "SRR1044416", "SRR1044417", "SRR1044418",
        "SRR1044419", "SRR1044420", "SRR1044421"
    ],
    "MYOCARDIUM_51nt/02_FETAL": [
        "SRR1044513", "SRR1044515"
    ],
    "MYOCARDIUM_51nt/03_DCM_ADVANCED": [
        "SRR1044452", "SRR1044454", "SRR1044457", "SRR1044458",
        "SRR1044461", "SRR1044469", "SRR1044474", "SRR1044475",
        "SRR1044476", "SRR1044477", "SRR1044478", "SRR1044479",
        "SRR1044482", "SRR1044485", "SRR1044487", "SRR1044489",
        "SRR1044494", "SRR1044495", "SRR1044496", "SRR1044500",
        "SRR1044502"
    ],
    "MYOCARDIUM_51nt/04_DCM_EXPLANT": [
        "SRR1044451", "SRR1044453", "SRR1044460", "SRR1044468",
        "SRR1044488", "SRR1044499", "SRR1044501"
    ],
    "MYOCARDIUM_51nt/05_ICM_ADVANCED": [
        "SRR1044527", "SRR1044524", "SRR1044529", "SRR1044535",
        "SRR1044536", "SRR1044544", "SRR1044549", "SRR1044551",
        "SRR1044555", "SRR1044560", "SRR1044566", "SRR1044577",
        "SRR1044578"
    ],
    "MYOCARDIUM_51nt/06_ICM_EXPLANT": [
        "SRR1044523", "SRR1044528", "SRR1044534", "SRR1044550",
        "SRR1044565", "SRR1044576"
    ],

    # MYOCARDIUM_36nt
    "MYOCARDIUM_36nt/01_NF": [
        "SRR1044422", "SRR1044423", "SRR1044424"
    ],
    "MYOCARDIUM_36nt/02_DCM_ADVANCED": [
        "SRR1044456"
    ],
    "MYOCARDIUM_36nt/03_DCM_EXPLANT": [
        "SRR1044455"
    ],
    "MYOCARDIUM_36nt/04_ICM_ADVANCED": [
        "SRR1044522"
    ],
    "MYOCARDIUM_36nt/05_ICM_EXPLANT": [
        "SRR1044521"
    ],

    # PLASMA_51nt
    "PLASMA_51nt/01_NF": [
        "SRR1044425", "SRR1044427", "SRR1044429", "SRR1044431",
        "SRR1044433", "SRR1044434", "SRR1044435", "SRR1044436",
        "SRR1044437", "SRR1044438", "SRR1044439", "SRR1044440",
        "SRR1044441", "SRR1044442", "SRR1044443", "SRR1044444",
        "SRR1044445"
    ],
    "PLASMA_51nt/02_DCM_ADVANCED": [
        "SRR1044464", "SRR1044465", "SRR1044471", "SRR1044481",
        "SRR1044484", "SRR1044486", "SRR1044491", "SRR1044498",
        "SRR1044505"
    ],
    "PLASMA_51nt/03_DCM_3M_LVAD": [
        "SRR1044462", "SRR1044483", "SRR1044503"
    ],
    "PLASMA_51nt/04_DCM_6M_LVAD": [
        "SRR1044480", "SRR1044497", "SRR1044504"
    ],
    "PLASMA_51nt/05_DCM_EXPLANT": [
        "SRR1044463", "SRR1044470", "SRR1044490"
    ],
    "PLASMA_51nt/06_ICM_ADVANCED": [
        "SRR1044526", "SRR1044531", "SRR1044539", "SRR1044540",
        "SRR1044546", "SRR1044548", "SRR1044554", "SRR1044557",
        "SRR1044559", "SRR1044562", "SRR1044564", "SRR1044568",
        "SRR1044569", "SRR1044573", "SRR1044580", "SRR1044582",
        "SRR1044584", "SRR1044586"
    ],
    "PLASMA_51nt/07_ICM_3M_LVAD": [
        "SRR1044537", "SRR1044547", "SRR1044552", "SRR1044563",
        "SRR1044579", "SRR1044581", "SRR1044585"
    ],
    "PLASMA_51nt/08_ICM_6M_LVAD": [
        "SRR1044525", "SRR1044545", "SRR1044553", "SRR1044556",
        "SRR1044558", "SRR1044561", "SRR1044583"
    ],
    "PLASMA_51nt/09_ICM_EXPLANT": [
        "SRR1044530", "SRR1044538", "SRR1044567", "SRR1044572"
    ],
    "PLASMA_51nt/10_STABLE_HF": [
        "SRR1044506", "SRR1044507", "SRR1044508", "SRR1044509",
        "SRR1044587", "SRR1044588", "SRR1044589", "SRR1044590",
        "SRR1044591", "SRR1044592", "SRR1044593", "SRR1044594",
        "SRR1044595", "SRR1044596"
    ],

    # SERUM_51nt
    "SERUM_51nt/01_NF": [
        "SRR1044426", "SRR1044428", "SRR1044430", "SRR1044432"
    ],
    "SERUM_51nt/02_DCM_ADVANCED": [
        "SRR1044467", "SRR1044473", "SRR1044493"
    ],
    "SERUM_51nt/03_DCM_EXPLANT": [
        "SRR1044466", "SRR1044472", "SRR1044492"
    ],
    "SERUM_51nt/04_ICM_ADVANCED": [
        "SRR1044533", "SRR1044543", "SRR1044571", "SRR1044575"
    ],
    "SERUM_51nt/05_ICM_EXPLANT": [
        "SRR1044532", "SRR1044542", "SRR1044570", "SRR1044574"
    ],

    # OTHER_TISSUES
    "OTHER_TISSUES/ADIPOSE": ["SRR1044446"],
    "OTHER_TISSUES/BRAIN": ["SRR1044447"],
    "OTHER_TISSUES/LIVER": ["SRR1044448"],
    "OTHER_TISSUES/SKIN": ["SRR1044449"],
    "OTHER_TISSUES/SPLEEN": ["SRR1044450"],
    "OTHER_TISSUES/SKELETAL_MUSCLE": [
        "SRR1044514", "SRR1044459", "SRR1044541", "SRR1044518"
    ],
    "OTHER_TISSUES/HUVEC": ["SRR1044519", "SRR1044520"],
    "OTHER_TISSUES/PBMC": ["SRR1044597"],
    "OTHER_TISSUES/RBC": ["SRR1044598", "SRR1044599"]
}

EXCLUDED_SAMPLES = {
    "SRR1044510": "76nt fetal cardiac",
    "SRR1044511": "76nt fetal cardiac",
    "SRR1044512": "76nt fetal cardiac",
    "SRR1044516": "36nt fetal cardiac",
    "SRR1044517": "76nt fetal cardiac"
}

CORRECTIONS_LOG = {
    "SRR1044528": "Moved from ICM_ADVANCED to ICM_EXPLANT only",
    "SRR1044530": "Moved from myocardium to PLASMA_09_ICM_EXPLANT only",
    "SRR1044574": "Moved from myocardium to SERUM_05_ICM_EXPLANT only",
    "SRR1044432": "Moved from plasma to SERUM_01_NF only",
    "SRR1044521": "New folder MYOCARDIUM_36nt/05_ICM_EXPLANT created"
}


def check_file_exists(folder_path, srr_id):
    """Check if SRR file exists with any valid extension."""
    extensions = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
    for ext in extensions:
        fpath = folder_path / f"{srr_id}{ext}"
        if fpath.exists():
            return fpath, fpath.stat().st_size
    return None, 0


def main():
    print("=" * 70)
    print("GSE53080 FOLDER STRUCTURE VERIFICATION")
    print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    base = Path(BASE_PATH)

    if not base.exists():
        print(f"\nERROR: Base path does not exist: {base}")
        print("Please update BASE_PATH in the script to your actual GSE53080 folder.")
        sys.exit(1)

    print(f"\nBase path: {base}")
    print(f"Total folders to check: {len(EXPECTED_STRUCTURE)}")
    print(f"Total expected files: {sum(len(v) for v in EXPECTED_STRUCTURE.values())}")

    # Statistics
    stats = {
        "folders_ok": 0,
        "folders_missing": 0,
        "folders_incomplete": 0,
        "files_expected": 0,
        "files_found": 0,
        "files_missing": [],
        "files_unexpected": [],
        "total_size_bytes": 0
    }

    # Check each folder
    for folder_rel, expected_srrs in EXPECTED_STRUCTURE.items():
        folder_path = base / folder_rel
        stats["files_expected"] += len(expected_srrs)

        print(f"\n{'─' * 60}")
        print(f"📁 {folder_rel}")
        print(f"   Expected: {len(expected_srrs)} files")

        if not folder_path.exists():
            stats["folders_missing"] += 1
            print(f"   ❌ FOLDER MISSING!")
            for srr in expected_srrs:
                stats["files_missing"].append((folder_rel, srr, "folder_missing"))
            continue

        # Check expected files
        found_count = 0
        for srr in expected_srrs:
            fpath, size = check_file_exists(folder_path, srr)
            if fpath:
                found_count += 1
                stats["files_found"] += 1
                stats["total_size_bytes"] += size
                size_mb = size / (1024**2)
                print(f"   ✅ {srr}.fastq.gz ({size_mb:.1f} MB)")
            else:
                stats["files_missing"].append((folder_rel, srr, "file_not_found"))
                print(f"   ❌ {srr}.fastq.gz MISSING")

        # Check for unexpected files
        if folder_path.exists():
            for item in folder_path.iterdir():
                if item.is_file():
                    srr_from_file = item.name.split('.')[0]
                    if srr_from_file not in expected_srrs:
                        stats["files_unexpected"].append((folder_rel, item.name))
                        print(f"   ⚠️  UNEXPECTED: {item.name}")

        # Folder status
        if found_count == len(expected_srrs):
            stats["folders_ok"] += 1
            print(f"   📊 Status: COMPLETE ({found_count}/{len(expected_srrs)})")
        else:
            stats["folders_incomplete"] += 1
            print(f"   📊 Status: INCOMPLETE ({found_count}/{len(expected_srrs)})")

    # Check excluded samples not present
    print(f"\n{'=' * 60}")
    print("EXCLUDED SAMPLES CHECK (should NOT be in analysis folders)")
    print(f"{'=' * 60}")
    excluded_found = []
    for srr, reason in EXCLUDED_SAMPLES.items():
        found_in = []
        for folder_rel in EXPECTED_STRUCTURE.keys():
            folder_path = base / folder_rel
            if folder_path.exists():
                fpath, _ = check_file_exists(folder_path, srr)
                if fpath:
                    found_in.append(folder_rel)
        if found_in:
            excluded_found.append((srr, reason, found_in))
            print(f"⚠️  {srr} ({reason}) FOUND in: {found_in}")
        else:
            print(f"✅ {srr} ({reason}) — correctly excluded")

    # Verify corrections
    print(f"\n{'=' * 60}")
    print("CORRECTIONS VERIFICATION")
    print(f"{'=' * 60}")
    for srr, correction in CORRECTIONS_LOG.items():
        locations = []
        for folder_rel, srrs in EXPECTED_STRUCTURE.items():
            if srr in srrs:
                folder_path = base / folder_rel
                if folder_path.exists():
                    fpath, _ = check_file_exists(folder_path, srr)
                    if fpath:
                        locations.append(folder_rel)

        if len(locations) == 1:
            print(f"✅ {srr}: {correction}")
            print(f"   Found in exactly 1 folder: {locations[0]}")
        elif len(locations) == 0:
            print(f"❌ {srr}: {correction}")
            print(f"   NOT FOUND anywhere!")
        else:
            print(f"❌ {srr}: {correction}")
            print(f"   DUPLICATE in: {locations} — should be in exactly 1 folder!")

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY REPORT")
    print(f"{'=' * 70}")
    print(f"Folders OK:        {stats['folders_ok']}/{len(EXPECTED_STRUCTURE)}")
    print(f"Folders Missing:   {stats['folders_missing']}")
    print(f"Folders Incomplete: {stats['folders_incomplete']}")
    print(f"Files Expected:    {stats['files_expected']}")
    print(f"Files Found:       {stats['files_found']}")
    print(f"Files Missing:     {len(stats['files_missing'])}")
    print(f"Unexpected Files:  {len(stats['files_unexpected'])}")
    print(f"Total Size:        {stats['total_size_bytes'] / (1024**3):.2f} GB")
    print(f"Excluded Found:    {len(excluded_found)} (should be 0)")

    if stats['files_missing']:
        print(f"\n{'─' * 60}")
        print("MISSING FILES LIST:")
        for folder, srr, reason in stats['files_missing']:
            print(f"  {folder}/{srr}.fastq.gz ({reason})")

    if stats['files_unexpected']:
        print(f"\n{'─' * 60}")
        print("UNEXPECTED FILES LIST:")
        for folder, fname in stats['files_unexpected']:
            print(f"  {folder}/{fname}")

    # Final verdict
    print(f"\n{'=' * 70}")
    if (stats['folders_ok'] == len(EXPECTED_STRUCTURE) and 
        len(stats['files_missing']) == 0 and 
        len(excluded_found) == 0):
        print("✅ ALL CHECKS PASSED — Ready for analysis!")
    else:
        print("❌ ISSUES FOUND — Please fix before proceeding!")
    print(f"{'=' * 70}")

    # Save report to file
    report_file = base / "verification_report.txt"
    with open(report_file, 'w') as f:
        f.write(f"GSE53080 Verification Report\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Base: {base}\n\n")
        f.write(f"Folders OK: {stats['folders_ok']}/{len(EXPECTED_STRUCTURE)}\n")
        f.write(f"Files Found: {stats['files_found']}/{stats['files_expected']}\n")
        f.write(f"Total Size: {stats['total_size_bytes'] / (1024**3):.2f} GB\n")
        if stats['files_missing']:
            f.write(f"\nMissing Files:\n")
            for folder, srr, reason in stats['files_missing']:
                f.write(f"  {folder}/{srr}.fastq.gz\n")

    print(f"\nReport saved to: {report_file}")


if __name__ == "__main__":
    main()
