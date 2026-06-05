# GSE53080_pipeline

Pipeline for analyzing small RNA-seq data from GSE53080 (Akat et al. 2014, PNAS).

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

