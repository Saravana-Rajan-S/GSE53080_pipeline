# GSE53080_pipeline

Pipeline for analyzing small RNA-seq data from GSE53080 (Akat et al. 2014, PNAS).

**Dataset:** Akat et al. 2014, PNAS 111(30):11151-11156  
**GEO Accession:** GSE53080 | **SRA Study:** SRP033566  
**Total samples:** 185 (metadata) → 180 analyzed (5 excluded for batch effect)

---

## 📁 Folder Structure

The GSE53080 dataset is organized into 35 folders based on tissue type, read length, and disease status. This structure was derived from metadata cross-checking against the SraRunTable (PRJNA230811) and includes 5 critical corrections.

### Organization Logic

```
GSE53080/
├── MYOCARDIUM_51nt/          # Heart tissue, 51nt reads (main analysis)
│   ├── 01_NF/                # Non-failing controls (7 samples)
│   ├── 02_FETAL/             # Fetal heart (2 samples)
│   ├── 03_DCM_ADVANCED/      # Dilated cardiomyopathy, advanced (21)
│   ├── 04_DCM_EXPLANT/       # DCM, post-LVAD explant (7)
│   ├── 05_ICM_ADVANCED/      # Ischemic cardiomyopathy, advanced (13)
│   └── 06_ICM_EXPLANT/       # ICM, post-LVAD explant (6)
├── MYOCARDIUM_36nt/          # Heart tissue, 36nt reads (subset)
│   ├── 01_NF/                # Non-failing (3)
│   ├── 02_DCM_ADVANCED/      # DCM advanced (1)
│   ├── 03_DCM_EXPLANT/       # DCM explant (1)
│   ├── 04_ICM_ADVANCED/      # ICM advanced (1)
│   └── 05_ICM_EXPLANT/       # ICM explant (1)
├── PLASMA_51nt/              # Blood plasma, 51nt reads
│   ├── 01_NF/                # Non-failing (17)
│   ├── 02_DCM_ADVANCED/      # DCM advanced (9)
│   ├── 03_DCM_3M_LVAD/       # DCM, 3 months LVAD (3)
│   ├── 04_DCM_6M_LVAD/       # DCM, 6 months LVAD (3)
│   ├── 05_DCM_EXPLANT/       # DCM explant (3)
│   ├── 06_ICM_ADVANCED/      # ICM advanced (18)
│   ├── 07_ICM_3M_LVAD/       # ICM, 3 months LVAD (7)
│   ├── 08_ICM_6M_LVAD/       # ICM, 6 months LVAD (7)
│   ├── 09_ICM_EXPLANT/       # ICM explant (4)
│   └── 10_STABLE_HF/         # Stable heart failure (14)
├── SERUM_51nt/               # Blood serum, 51nt reads
│   ├── 01_NF/                # Non-failing (4)
│   ├── 02_DCM_ADVANCED/      # DCM advanced (3)
│   ├── 03_DCM_EXPLANT/       # DCM explant (3)
│   ├── 04_ICM_ADVANCED/      # ICM advanced (4)
│   └── 05_ICM_EXPLANT/       # ICM explant (4)
└── OTHER_TISSUES/            # Non-cardiac tissues (13 samples)
    ├── ADIPOSE, BRAIN, LIVER, SKIN, SPLEEN
    ├── SKELETAL_MUSCLE, HUVEC, PBMC, RBC
```

**Total:** 180 samples assigned | 5 excluded (fetal batch effect) | 185 in metadata

### Key Corrections Applied

| Sample | Original Assignment | Corrected Assignment | Reason |
|--------|--------------------|---------------------|--------|
| SRR1044528 | MYOCARDIUM_51nt/05_ICM_ADVANCED | MYOCARDIUM_51nt/06_ICM_EXPLANT | `lvad=Explantation` in metadata |
| SRR1044530 | MYOCARDIUM_51nt/06_ICM_EXPLANT | PLASMA_51nt/09_ICM_EXPLANT | `tissue=Plasma` in metadata |
| SRR1044574 | MYOCARDIUM_51nt/06_ICM_EXPLANT | SERUM_51nt/05_ICM_EXPLANT | `tissue=Serum` in metadata |
| SRR1044432 | PLASMA_51nt/01_NF | SERUM_51nt/01_NF | `tissue=Serum` in metadata |
| SRR1044521 | MYOCARDIUM_36nt/04_ICM_ADVANCED | MYOCARDIUM_36nt/05_ICM_EXPLANT | `lvad=Explantation` in metadata |

> **Full details:** See [`docs/folder_assignments.md`](docs/folder_assignments.md)

---

## 🔬 Barcode Verification (COMPLETE)

### Adapter Structure

```
[5nt sample-specific barcode] + TCGTATGCCGTCTTCTGCTTG
         ↑                        ↑
      5nt barcode           4nt invariant (TCGT)
      (20 unique)           + 17nt Illumina adapter
```

- **20 adapters (Ad01–Ad20)** with unique 5nt barcodes
- **Full adapter:** 26nt
- **9nt signature:** 5nt barcode + TCGT (used for verification)
- **Match policy:** EXACT MATCH, 0 mismatches
- **Search method:** SLIDING SEARCH through full read (no fixed position)

### Verification Results (180 files)

| Status | Count | % | Tissue Distribution |
|--------|-------|---|---------------------|
| **Clean** | 152 | 84.4% | All myocardium, other tissues, most plasma |
| **WARNING** | 28 | 15.6% | Plasma (16) + Serum (12) only |
| INFO | 0 | 0% | — |
| CRITICAL | 0 | 0% | — |

### Adapter Rate by Tissue

| Tissue | Avg Adapter Rate | Assessment |
|--------|-----------------|------------|
| Myocardium (51nt) | 99.4% | Excellent |
| Myocardium (36nt) | 96.7% | Excellent |
| Other tissues | 98.9% | Excellent |
| Plasma (51nt) | 91.1% | Good |
| Serum (51nt) | 71.6% | Lower (biologically expected) |

### Key Metrics

- **Barcode start position median:** 23 (1-based) → 22nt RNA insert, typical for mature miRNAs
- **Read lengths:** 168×51nt, 8×36nt, 2×40nt, 2×76nt
- **All barcodes cleanly assigned:** No mixed barcodes → demultiplexing was correct

**Files:**
- [`docs/barcode_verification_report.csv`](docs/barcode_verification_report.csv) — Full 180-sample report
- [`config/barcode_mapping.csv`](config/barcode_mapping.csv) — Clean SRR → Barcode mapping for trimming
- [`src/preprocessing/verify_barcodes.py`](src/preprocessing/verify_barcodes.py) — Verification script

---

## 🧪 28 WARNING Files — 5nt vs 9nt Diagnostic (COMPLETE)

### The Problem

28 files showed "low adapter rate" in the initial 9nt scan:
- **9nt rate:** 49.95% – 77.07%
- All are **plasma** or **serum** (no myocardium or other tissues)

**Hypothesis:** The 5nt barcode is present at 100%, but the TCGT suffix is degraded/mutated in ~30–50% of reads. This is biologically expected in plasma/serum due to RNase activity.

### Diagnostic Method (v4)

1. **Single-pass sliding scan** — 9nt checked first; if not found, 5nt counted with `seq.count()`
2. **Double-hit detection** — counts reads where the 5nt barcode appears **2+ times** (once real adapter, once random noise inside RNA)
3. **Threshold:** ≥10% double-hit rate = significant false positive contamination

### Diagnostic Results

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **5nt-total %** | **100.0% in ALL 28 files** | Every read contains the assigned barcode |
| **9nt %** | 49.95% – 77.07% | True adapter rate (barcode + intact TCGT) |
| **Gap (5nt − 9nt)** | 22.06% – 50.05% | Reads with degraded TCGT suffix |
| **Double-hit %** | **0.02% – 2.57%** | Random false positive contamination |
| **False positive threshold** | <10% | **All 28 files pass** |

### Conclusion

> **TCGT ERROR confirmed in all 28 WARNING files.** The 5nt barcode is present and correctly placed in 100% of reads. The ~30–50% gap between 5nt and 9nt represents genuine adapter molecules where the TCGT suffix is degraded or mutated — **not** barcode misassignment or random noise.

**Why serum is worse than plasma:** Serum is collected after blood clotting, which releases RNases that degrade RNA and damage the adapter sequence. Plasma (anti-coagulated) avoids this. This is biologically expected and validates data quality.

### Trimming Strategy

**All 180 files are suitable for trimming at the 5nt barcode position.**

- Trim point: right before the 5nt barcode (wherever it occurs in the read)
- Tool: `cutadapt` with the full 26nt adapter sequence
- The 5nt barcode is the correct trim anchor for all scenarios

**Files:**
- [`docs/warning_files_5nt_vs_9nt_diagnostic.csv`](docs/warning_files_5nt_vs_9nt_diagnostic.csv) — Full 28-file diagnostic report
- [`src/preprocessing/verify_28_warning_files.py`](src/preprocessing/verify_28_warning_files.py) — Diagnostic script (v4)

---

## 🛠️ Tools & References

### Installed Tools
- `cutadapt` — adapter trimming
- `bowtie2` — alignment
- `samtools` — BAM manipulation
- `fastqc` — quality control
- `multiqc` — report aggregation
- `bedtools` — genomic operations
- `featureCounts` (subread) — read counting
- R + Bioconductor (DESeq2, edgeR, ggplot2)
- Python (pandas, numpy, pysam)

### Reference Databases
Located at `/home/saravana/Downloads/smallRNA_references/`:
- `01_mirbase/` — hairpin.fa, mature.fa
- `02_yrna/` — YRNA_reference.fa + .bt2 indices
- `03_gtrnadb/` — hg38-mature-tRNAs.fa
- `04_snorna/` — human_snoRNAs.fa
- `05_pirna/` — hsa.gold.fa, hsa.v3.0.fa, Cardiovascular.txt
- `06_hg38/` — 6 .bt2 index files
- `07_combined/` — smallRNA_combined_dedup.fa + 6 .bt2 index files

---

## 📋 Pipeline Steps

| Step | Status | Description |
|------|--------|-------------|
| 1. Folder structure | ✅ Complete | 180 samples organized, 5 corrections applied |
| 2. Barcode verification | ✅ Complete | 180 files scanned, 28 flagged for TCGT analysis |
| 3. 5nt vs 9nt diagnostic | ✅ Complete | All 28 files confirmed TCGT error, ready for trimming |
| 4. FastQC (pre-trim) | ⬜ Pending | Run on representative samples, aggregate with MultiQC |
| 5. Trimming (cutadapt) | ⬜ Pending | Trim at 5nt barcode position using full adapter |
| 6. Alignment (3-round) | ⬜ Pending | Bowtie2 to smallRNA index → hg38 → hairpin |
| 7. Quantification | ⬜ Pending | featureCounts / custom counting, build count matrix |
| 8. Normalization | ⬜ Pending | TPM, DESeq2 size factors |
| 9. Differential Expression | ⬜ Pending | DESeq2, edgeR, volcano plots |
| 10. tRF/yRF & isomiR | ⬜ Pending | Classification and extraction |
| 11. Machine Learning | ⬜ Pending | RF, XGBoost, LASSO, SHAP |

---

## 📝 Citation

Akat KM, Moore-McGriff D, Morozov P, Brown M, Gogakos T, Correa da Rosa J, et al. Comparative RNA-sequencing analysis of myocardial and circulating small RNAs in human heart failure and their utility as biomarkers. *Proc Natl Acad Sci USA*. 2014;111(30):11151-11156.

---

## 👤 Author

**Saravana Rajan S**  
GitHub: [Saravana-Rajan-S](https://github.com/Saravana-Rajan-S)
