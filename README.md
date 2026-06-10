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

### Trimming Strategy Applied

**All 180 files were trimmed at the 5nt barcode position using a dual-pass cutadapt strategy.**

- Pass 1: 9nt exact match (barcode + TCGT), 0 mismatches
- Pass 2: 5nt barcode rescue on Pass 1 leftovers
- Final output: `trimmed/5nt_matched/` — 99.89% read retention across all 180 samples

See the [Trimming section](#️-trimming--dual-pass-strategy-complete) for full results.

**Files:**
- [`docs/warning_files_5nt_vs_9nt_diagnostic.csv`](docs/warning_files_5nt_vs_9nt_diagnostic.csv) — Full 28-file diagnostic report
- [`src/preprocessing/verify_28_warning_files.py`](src/preprocessing/verify_28_warning_files.py) — Diagnostic script (v4)

---

## 🔍 FastQC / MultiQC — Pre-Trimming QC (COMPLETE)

### Overview

FastQC was run on all 180 FASTQ files with batch segregation by sequencing platform, then MultiQC aggregated each batch separately. **Zero samples were excluded.**

| Metric | Value |
|--------|-------|
| Files analysed | 180 |
| Tools | FastQC + MultiQC |
| Samples excluded | 0 (5 already excluded at folder structure step) |
| Batches | 3 (HiSeq2000 non-serum, HiSeq2000 serum, GA IIx) |

### Sequencing Platform Batches

Samples were segregated by sequencer from SraRunTable metadata before MultiQC aggregation — mixing platforms in a single MultiQC report obscures platform-specific quality artefacts.

| Batch | Platform | Samples | MultiQC Report |
|-------|----------|---------|----------------|
| A1 | Illumina HiSeq 2000 — non-serum | 132 | `qc/multiqc/batch_hiseq/` |
| A2 | Illumina HiSeq 2000 — serum | 18 | `qc/multiqc/batch_hiseq_serum/` |
| B | Illumina GA IIx | 30 | `qc/multiqc/batch_gaiix/` |

FastQC outputs are stored per-batch:
```
qc/
├── fastqc/
│   ├── hiseq2000/
│   │   ├── non_serum/   ← 132 FastQC reports
│   │   └── serum/       ← 18 FastQC reports
│   └── gaiix/
│       └── all/         ← 30 FastQC reports
└── multiqc/
    ├── batch_hiseq/         ← HiSeq2000 non-serum report
    ├── batch_hiseq_serum/   ← HiSeq2000 serum report
    └── batch_gaiix/         ← GA IIx report
```

### Read Length Distribution

| Read length | Files | Tissue |
|-------------|-------|--------|
| 51 nt | 168 | Myocardium, Plasma, Serum, most Other |
| 36 nt | 8 | Myocardium 36nt group |
| 40 nt | 2 | Other tissues |
| **76 nt** | **2** | **Skeletal muscle (SRR1044514, SRR1044518) — GA IIx batch** |

### SRR1044514 & SRR1044518 — FastQC FAIL Flags Explained

The only two 76 nt files in the dataset (GA IIx batch, skeletal muscle). FastQC flagged multiple FAIL modules — all are artefacts, not data quality problems.

**Root cause:** The sequencer ran 76 cycles on ~22 nt inserts. After insert (22 nt) + adapter (26 nt) = 48 nt of meaningful content, the remaining ~28 positions read random flanking sequence — **adapter read-through**.

| FastQC FAIL | Cause | Action |
|-------------|-------|--------|
| Per base sequence quality | Quality drops at position ~48 — read-through zone | Trimming removes it |
| Per tile sequence quality | Tile-level artefact from same read-through | Trimming removes it |
| Per base sequence content | Chaotic %ATGC at positions 48–76 | Trimming removes it |
| Sequence Duplication Levels | ~2,500 unique miRNA sequences — duplication expected for small RNA | Ignore (module irrelevant for small RNA) |
| Overrepresented sequences | Real miRNA reads carrying barcode + adapter tail | Expected and correct |
| Adapter Content | Adapter genuinely present at position ~22 | Cutadapt handles it |

**Verdict:** Both files are clean skeletal muscle small RNA libraries. FAILs are 100% artefacts of a 76 nt run on 22 nt inserts. Trimming parameters are identical to all other files.

### Low-Read Files (<1M Reads)

23 files fell below 1M reads — all investigated and retained:

| Tissue type | Count | Explanation |
|-------------|-------|-------------|
| Plasma | 19 | Circulating miRNA packaged in exosomes/AGO2; 100–1000× less starting material than tissue |
| Serum | 3 | Additional RNase exposure during 30–60 min clotting; no EDTA anticoagulant |
| Liver (SRR1044448) | 1 | 989K reads — borderline, Clean status, technical variability |

Note: EDTA (present in plasma collection tubes) chelates Mg²⁺, inhibiting RNases — this is why plasma yields are generally better than serum within the same disease group.

### Trimming Readiness

| Check | Status |
|-------|--------|
| FastQC on all 180 files | ✅ |
| MultiQC — 3 batch reports | ✅ |
| SRR1044514/18 FAILs explained | ✅ Artefact, not defect |
| Low-read files confirmed biological | ✅ |
| Samples to exclude | ✅ None |

**Files:**
- [`src/preprocessing/run_fastqc_multiqc.py`](src/preprocessing/run_fastqc_multiqc.py) — FastQC + MultiQC runner with platform batch segregation
- [`qc/`](qc/) — FastQC HTML reports and MultiQC summaries (3 batch reports)



### Strategy

A two-pass cutadapt pipeline was designed to maximise read retention, specifically to rescue reads from plasma and serum samples where the TCGT suffix is degraded.

```
Pass 1 — 9nt exact match (Strategy A)
  cutadapt -a [barcode+TCGT] -O 9 -e 0 -m 15
  → Trimmed reads  → 9nt_match/
  → Untrimmed reads (TCGT-degraded) → rescued_input/

Pass 2 — 5nt rescue on leftovers (Strategy B)
  cutadapt -a [barcode] -O 5 -e 0 -m 15 --discard-untrimmed
  → Rescued reads  → 5nt_rescued/

Pass 3 — Merge
  cat 9nt_match/ + 5nt_rescued/ → 5nt_matched/   (final output)
```

> **Key design decision:** `-M 45` (maximum length filter) was intentionally removed from both passes. It was silently discarding reads before they could reach `--untrimmed-output`, causing reads to vanish without a rescue path. Length filtering is not applied at the trimming stage.

### Output Folders

| Folder | Contents | Files |
|--------|----------|-------|
| `trimmed/9nt_match/` | Pass 1 trimmed reads (9nt exact) | 180 |
| `trimmed/rescued_input/` | Untrimmed leftovers from Pass 1 | 180 |
| `trimmed/5nt_rescued/` | Pass 2 rescued reads (5nt match) | 180 |
| `trimmed/5nt_matched/` | Final merged output | 180 |

### Results

| Metric | Value |
|--------|-------|
| Total original reads | 917,831,939 |
| 9nt matched | 866,547,789 (94.4%) |
| 5nt rescued | 50,302,602 (5.5%) |
| **Final retained** | **916,850,391 (99.89%)** |
| Unrecoverable loss | 981,548 (0.11%) |
| Avg rescue % per sample | 7.33% |

### Rescue Rate by Sample Type

| Sample type | Avg rescue % | Notes |
|-------------|-------------|-------|
| Myocardium 51nt | ~0.5% | High quality, minimal degradation |
| Myocardium 36nt | ~3–5% | Shorter reads, more partial adapters |
| Plasma 51nt | 2–38% | Cell-free miRNA, moderate degradation |
| Serum 51nt | 22–34% | Post-clotting RNase activity, most degraded |

### Post-Trimming Verification

All 180 samples cross-validated against the pre-trimming barcode report:
- Barcode, Adapter ID, Tissue group, Total reads, Adapter %, Status — **all 180/180 match**
- RNA insert length median post-trim = pre-trim median — **180/180 exact**
- Residual 5nt barcode signal: 1.89% — confirmed as internal biological false positives (random 5-mer occurrence), not trimming failures

**Files:**
- [`src/preprocessing/run_dual_trimming.py`](src/preprocessing/run_dual_trimming.py) — Dual-pass trimming pipeline
- [`src/preprocessing/verify_post_trimming.py`](src/preprocessing/verify_post_trimming.py) — Post-trimming verification script
- [`docs/trimming_comparison_summary.csv`](docs/trimming_comparison_summary.csv) — Per-sample comparison (180 rows)
- [`docs/post_trimming_verification.csv`](docs/post_trimming_verification.csv) — Per-file verification (720 rows)



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
| 4. FastQC (pre-trim) | ✅ Complete | All 180 files; MultiQC aggregated; 0 samples excluded |
| 5. Trimming (cutadapt) | ✅ Complete | Dual-pass strategy: 9nt exact + 5nt rescue; 99.89% retention |
| 6. Post-trim verification | ✅ Complete | 180/180 verified; 50.3M reads rescued (5.48%) |
| 7. Alignment (3-round) | ⬜ Pending | Bowtie2 to smallRNA index → hg38 → hairpin |
| 8. Quantification | ⬜ Pending | featureCounts / custom counting, build count matrix |
| 9. Normalization | ⬜ Pending | TPM, DESeq2 size factors |
| 10. Differential Expression | ⬜ Pending | DESeq2, edgeR, volcano plots |
| 11. tRF/yRF & isomiR | ⬜ Pending | Classification and extraction |
| 12. Machine Learning | ⬜ Pending | RF, XGBoost, LASSO, SHAP |

---

## 📝 Citation

Akat KM, Moore-McGriff D, Morozov P, Brown M, Gogakos T, Correa da Rosa J, et al. Comparative RNA-sequencing analysis of myocardial and circulating small RNAs in human heart failure and their utility as biomarkers. *Proc Natl Acad Sci USA*. 2014;111(30):11151-11156.

---

## 👤 Author

**Saravana Rajan S**  
GitHub: [Saravana-Rajan-S](https://github.com/Saravana-Rajan-S)
