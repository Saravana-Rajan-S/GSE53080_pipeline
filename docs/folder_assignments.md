# GSE53080 — Final Verified Folder Assignments

**Dataset:** Akat et al. 2014, PNAS 111(30):11151-11156  
**GEO Accession:** GSE53080 | **SRA Study:** SRP033566  
**Generated from:** SraRunTable metadata (PRJNA230811) — FULLY CROSS-CHECKED  
**Date:** 2026-05-29

---

## Corrections Applied (based on metadata verification)

1. **SRR1044528** removed from `MYOCARDIUM_51nt/05_ICM_ADVANCED` → moved to `MYOCARDIUM_51nt/06_ICM_EXPLANT` (metadata: `lvad=Explantation`)
2. **SRR1044530** removed from `MYOCARDIUM_51nt/06_ICM_EXPLANT` → moved to `PLASMA_51nt/09_ICM_EXPLANT` (metadata: `tissue=Plasma`)
3. **SRR1044574** removed from `MYOCARDIUM_51nt/06_ICM_EXPLANT` → moved to `SERUM_51nt/05_ICM_EXPLANT` (metadata: `tissue=Serum`)
4. **SRR1044432** removed from `PLASMA_51nt/01_NF` → moved to `SERUM_51nt/01_NF` (metadata: `tissue=Serum`)
5. **MYOCARDIUM_36nt/04_ICM_ADVANCED** split into `04_ICM_ADVANCED` + `05_ICM_EXPLANT`  
   (SRR1044521 has `lvad=Explantation`; mirrors 51nt folder structure for pipeline compatibility)

## Excluded from Analysis (read-length batch incompatibility)

| SRR ID | Reason |
|--------|--------|
| SRR1044510 | 76nt fetal cardiac muscle |
| SRR1044511 | 76nt fetal cardiac muscle |
| SRR1044512 | 76nt fetal cardiac muscle |
| SRR1044516 | 36nt fetal cardiac muscle |
| SRR1044517 | 76nt fetal cardiac muscle |

---

## Folder Structure


### `MYOCARDIUM_51nt/01_NF` — 7 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044415 | age=60.0 | gender=male | lvad=No LVAD |
| SRR1044416 | age=60.0 | gender=male | lvad=No LVAD |
| SRR1044417 | age=56.0 | gender=male | lvad=No LVAD |
| SRR1044418 | age=56.0 | gender=male | lvad=No LVAD |
| SRR1044419 | age=40.0 | gender=male | lvad=No LVAD |
| SRR1044420 | age=2.0 | gender=male | lvad=No LVAD |
| SRR1044421 | age=80.0 | gender=male | lvad=No LVAD |

### `MYOCARDIUM_51nt/02_FETAL` — 2 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044513 | age=0.0 | gender=female | lvad=No LVAD |
| SRR1044515 | age=0.0 | lvad=No LVAD |

### `MYOCARDIUM_51nt/03_DCM_ADVANCED` — 21 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044452 | age=72.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044454 | age=44.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044457 | age=44.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044458 | age=44.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044461 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044469 | age=52.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044474 | age=46.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044475 | age=62.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044476 | age=48.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044477 | age=65.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044478 | age=57.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044479 | age=78.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044482 | age=51.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044485 | age=74.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044487 | age=72.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044489 | age=66.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044494 | age=40.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044495 | age=70.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044496 | age=33.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044500 | age=45.0 | gender=female | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044502 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |

### `MYOCARDIUM_51nt/04_DCM_EXPLANT` — 7 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044451 | age=72.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044453 | age=44.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044460 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044468 | age=52.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044488 | age=66.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044499 | age=45.0 | gender=female | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044501 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |

### `MYOCARDIUM_51nt/05_ICM_ADVANCED` — 13 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044527 | age=78.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044524 | age=62.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044529 | age=70.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044535 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044536 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044544 | age=54.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044549 | age=73.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044551 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044555 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044560 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044566 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044577 | age=68.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044578 | age=69.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |

### `MYOCARDIUM_51nt/06_ICM_EXPLANT` — 6 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044523 | age=62.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044528 | age=70.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044534 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044550 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044565 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044576 | age=68.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |

### `MYOCARDIUM_36nt/01_NF` — 3 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044422 | lvad=No LVAD |
| SRR1044423 | lvad=No LVAD |
| SRR1044424 | lvad=No LVAD |

### `MYOCARDIUM_36nt/02_DCM_ADVANCED` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044456 | age=58.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |

### `MYOCARDIUM_36nt/03_DCM_EXPLANT` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044455 | age=58.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |

### `MYOCARDIUM_36nt/04_ICM_ADVANCED` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044522 | age=51.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |

### `MYOCARDIUM_36nt/05_ICM_EXPLANT` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044521 | age=51.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |

### `PLASMA_51nt/01_NF` — 17 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044425 | age=51.0 | gender=male | lvad=No LVAD |
| SRR1044427 | age=41.0 | gender=female | lvad=No LVAD |
| SRR1044429 | age=49.0 | gender=female | lvad=No LVAD |
| SRR1044431 | age=69.0 | gender=female | lvad=No LVAD |
| SRR1044433 | age=62.0 | gender=female | lvad=No LVAD |
| SRR1044434 | age=49.0 | gender=male | lvad=No LVAD |
| SRR1044435 | age=68.0 | gender=male | lvad=No LVAD |
| SRR1044436 | age=68.0 | gender=male | lvad=No LVAD |
| SRR1044437 | age=45.0 | gender=male | lvad=No LVAD |
| SRR1044438 | age=60.0 | gender=male | lvad=No LVAD |
| SRR1044439 | age=60.0 | gender=male | lvad=No LVAD |
| SRR1044440 | age=32.0 | gender=male | lvad=No LVAD |
| SRR1044441 | age=61.0 | gender=male | lvad=No LVAD |
| SRR1044442 | age=61.0 | gender=male | lvad=No LVAD |
| SRR1044443 | age=70.0 | gender=male | lvad=No LVAD |
| SRR1044444 | age=61.0 | gender=male | lvad=No LVAD |
| SRR1044445 | age=61.0 | gender=male | lvad=No LVAD |

### `PLASMA_51nt/02_DCM_ADVANCED` — 9 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044464 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044465 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044471 | age=52.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044481 | age=78.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044484 | age=51.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044486 | age=74.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044491 | age=66.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044498 | age=33.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044505 | age=49.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |

### `PLASMA_51nt/03_DCM_3M_LVAD` — 3 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044462 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044483 | age=51.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044503 | age=49.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=3 months LVAD support |

### `PLASMA_51nt/04_DCM_6M_LVAD` — 3 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044480 | age=78.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044497 | age=33.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044504 | age=49.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=6 months LVAD support |

### `PLASMA_51nt/05_DCM_EXPLANT` — 3 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044463 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044470 | age=52.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044490 | age=66.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |

### `PLASMA_51nt/06_ICM_ADVANCED` — 18 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044526 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044531 | age=70.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044539 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044540 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044546 | age=54.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044548 | age=73.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044554 | age=74.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044557 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044559 | age=68.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044562 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044564 | age=61.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044568 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044569 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044573 | age=70.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044580 | age=72.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044582 | age=75.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044584 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044586 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |

### `PLASMA_51nt/07_ICM_3M_LVAD` — 7 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044537 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044547 | age=73.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044552 | age=74.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044563 | age=61.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044579 | age=72.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044581 | age=75.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |
| SRR1044585 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=3 months LVAD support |

### `PLASMA_51nt/08_ICM_6M_LVAD` — 7 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044525 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044545 | age=54.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044553 | age=74.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044556 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044558 | age=68.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044561 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |
| SRR1044583 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=6 months LVAD support |

### `PLASMA_51nt/09_ICM_EXPLANT` — 4 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044530 | age=70.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044538 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044567 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044572 | age=70.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |

### `PLASMA_51nt/10_STABLE_HF` — 14 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044506 | age=49.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044507 | age=63.0 | gender=female | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044508 | age=63.0 | gender=female | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044509 | age=70.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044587 | age=62.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044588 | age=52.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044589 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044590 | age=63.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044591 | age=67.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044592 | age=61.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044593 | age=42.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044594 | age=65.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044595 | age=51.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044596 | age=71.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |

### `SERUM_51nt/01_NF` — 4 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044426 | age=51.0 | gender=male | lvad=No LVAD |
| SRR1044428 | age=41.0 | gender=female | lvad=No LVAD |
| SRR1044430 | age=49.0 | gender=female | lvad=No LVAD |
| SRR1044432 | age=69.0 | gender=female | lvad=No LVAD |

### `SERUM_51nt/02_DCM_ADVANCED` — 3 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044467 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044473 | age=52.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |
| SRR1044493 | age=66.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=No LVAD |

### `SERUM_51nt/03_DCM_EXPLANT` — 3 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044466 | age=61.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044472 | age=52.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |
| SRR1044492 | age=66.0 | gender=male | disease=Dilated Cardiomyopathy | lvad=Explantation |

### `SERUM_51nt/04_ICM_ADVANCED` — 4 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044533 | age=70.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044543 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044571 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |
| SRR1044575 | age=70.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=No LVAD |

### `SERUM_51nt/05_ICM_EXPLANT` — 4 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044532 | age=70.0 | gender=female | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044542 | age=66.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044570 | age=64.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |
| SRR1044574 | age=70.0 | gender=male | disease=Ischemic Cardiomyopathy | lvad=Explantation |

### `OTHER_TISSUES/ADIPOSE` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044446 | lvad=No LVAD |

### `OTHER_TISSUES/BRAIN` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044447 | lvad=No LVAD |

### `OTHER_TISSUES/LIVER` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044448 | lvad=No LVAD |

### `OTHER_TISSUES/SKIN` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044449 | lvad=No LVAD |

### `OTHER_TISSUES/SPLEEN` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044450 | lvad=No LVAD |

### `OTHER_TISSUES/SKELETAL_MUSCLE` — 4 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044514 | age=0.0 | gender=female | lvad=No LVAD |
| SRR1044459 | age=44.0 | gender=male | lvad=No LVAD |
| SRR1044541 | age=66.0 | gender=male | lvad=No LVAD |
| SRR1044518 | age=0.0 | lvad=No LVAD |

### `OTHER_TISSUES/HUVEC` — 2 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044519 | lvad=No LVAD |
| SRR1044520 | lvad=No LVAD |

### `OTHER_TISSUES/PBMC` — 1 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044597 | lvad=No LVAD |

### `OTHER_TISSUES/RBC` — 2 files

| SRR ID | Metadata |
|--------|----------|
| SRR1044598 | lvad=No LVAD |
| SRR1044599 | lvad=No LVAD |


---

## Summary

| Metric | Value |
|--------|-------|
| Total files assigned | 180 |
| Total in metadata | 185 |
| Excluded (batch effect) | 5 |
| Unassigned | 0 |
| Total folders | 35 |

---

*This document is the ground truth for the GSE53080 folder structure. Any future user needs this to reproduce the organization.*
