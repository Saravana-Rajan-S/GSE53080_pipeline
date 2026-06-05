#!/bin/bash
# ============================================================
# GSE53080 Quick Bash Verification
# Location in repo: src/utils/verify_gse53080.sh
#
# Run from repo root:
#   bash src/utils/verify_gse53080.sh /path/to/GSE53080
#
# Or with default path:
#   bash src/utils/verify_gse53080.sh
# ============================================================
# ============================================================

BASE="${1:-/home/saravana/Downloads/GSE53080}"

echo "=========================================="
echo "GSE53080 Quick Verification"
echo "Base: $BASE"
echo "=========================================="

# Check key folders exist
echo ""
echo "Checking folder structure..."
folders=(
    "MYOCARDIUM_51nt/01_NF" "MYOCARDIUM_51nt/02_FETAL" "MYOCARDIUM_51nt/03_DCM_ADVANCED"
    "MYOCARDIUM_51nt/04_DCM_EXPLANT" "MYOCARDIUM_51nt/05_ICM_ADVANCED" "MYOCARDIUM_51nt/06_ICM_EXPLANT"
    "MYOCARDIUM_36nt/01_NF" "MYOCARDIUM_36nt/02_DCM_ADVANCED" "MYOCARDIUM_36nt/03_DCM_EXPLANT"
    "MYOCARDIUM_36nt/04_ICM_ADVANCED" "MYOCARDIUM_36nt/05_ICM_EXPLANT"
    "PLASMA_51nt/01_NF" "PLASMA_51nt/02_DCM_ADVANCED" "PLASMA_51nt/03_DCM_3M_LVAD"
    "PLASMA_51nt/04_DCM_6M_LVAD" "PLASMA_51nt/05_DCM_EXPLANT" "PLASMA_51nt/06_ICM_ADVANCED"
    "PLASMA_51nt/07_ICM_3M_LVAD" "PLASMA_51nt/08_ICM_6M_LVAD" "PLASMA_51nt/09_ICM_EXPLANT"
    "PLASMA_51nt/10_STABLE_HF"
    "SERUM_51nt/01_NF" "SERUM_51nt/02_DCM_ADVANCED" "SERUM_51nt/03_DCM_EXPLANT"
    "SERUM_51nt/04_ICM_ADVANCED" "SERUM_51nt/05_ICM_EXPLANT"
    "OTHER_TISSUES/ADIPOSE" "OTHER_TISSUES/BRAIN" "OTHER_TISSUES/LIVER"
    "OTHER_TISSUES/SKIN" "OTHER_TISSUES/SPLEEN" "OTHER_TISSUES/SKELETAL_MUSCLE"
    "OTHER_TISSUES/HUVEC" "OTHER_TISSUES/PBMC" "OTHER_TISSUES/RBC"
)

missing_folders=0
for f in "${folders[@]}"; do
    if [ ! -d "$BASE/$f" ]; then
        echo "❌ MISSING FOLDER: $f"
        missing_folders=$((missing_folders + 1))
    fi
done

if [ $missing_folders -eq 0 ]; then
    echo "✅ All 35 folders present"
else
    echo "❌ $missing_folders folders missing"
fi

# Count files per folder
echo ""
echo "File counts per folder:"
for f in "${folders[@]}"; do
    if [ -d "$BASE/$f" ]; then
        count=$(ls -1 "$BASE/$f"/*.fastq.gz 2>/dev/null | wc -l)
        echo "  $f: $count files"
    fi
done

# Total fastq files
total=$(find "$BASE" -name "*.fastq.gz" | wc -l)
echo ""
echo "Total .fastq.gz files: $total (expected: 180)"

# Check for excluded samples (should NOT be present)
echo ""
echo "Checking excluded samples not present..."
excluded=("SRR1044510" "SRR1044511" "SRR1044512" "SRR1044516" "SRR1044517")
ex_found=0
for srr in "${excluded[@]}"; do
    found=$(find "$BASE" -name "$srr*" | wc -l)
    if [ $found -gt 0 ]; then
        echo "⚠️  $srr found (should be excluded)"
        ex_found=$((ex_found + 1))
    fi
done
if [ $ex_found -eq 0 ]; then
    echo "✅ All excluded samples correctly absent"
fi

# Check total size
echo ""
echo "Total dataset size:"
du -sh "$BASE"

echo ""
echo "=========================================="
if [ $missing_folders -eq 0 ] && [ $total -eq 180 ] && [ $ex_found -eq 0 ]; then
    echo "✅ ALL CHECKS PASSED"
else
    echo "❌ ISSUES FOUND - Review above"
fi
echo "=========================================="
