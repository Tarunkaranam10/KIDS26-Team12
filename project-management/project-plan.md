# Project Plan

## Goal

By the end of the event, Team 12 will determine whether technically harmonized DNA methylation can predict independently measured continuous reference HRDsum in unseen cancer types. The minimum result is an auditable TCGA cohort, frozen 450K/EPIC feature list, leakage-safe elastic-net baseline, at least one cancer-held-out evaluation, and a training-only-mean null comparator. A reproducible negative transfer result satisfies the scientific goal.

The strong result adds full adult multi-cancer transfer, a locked GBM/LGG evaluation, approved PBTP application, OOD/uncertainty reporting, interpretation, and an offline Shiny demonstration. CNV, fusion, component models, and augmentation remain stretch work.

## Tools

- R 4.5.0 or Python 3.10+ standard library for acquisition, checksums, sample matching, cohort QC, beta extraction, and engineering tests.
- R with `glmnet`, `data.table`, `jsonlite`, and `renv` for nested elastic-net training and frozen inference.
- Bioconductor/minfi and conumee2 only if the raw-intensity/CNV gate passes.
- Shiny as a presentation layer for approved precomputed results.
- Git and GitHub for version control; Slack for coordination.

## First Tasks

- [ ] Confirm the reported completed TCGA download path, expected size, checksum, and free storage; update `docs/14_VERIFICATION_STATUS.md` — Evan
- [ ] Install/activate R, run `tests/smoke_model.R`, and snapshot the tested environment — Evan, reviewed by Susanna
- [ ] Extract the eligible samples through `config/shared_autosomal_probes.txt`; rerun cohort QC and inspect per-cancer counts/missingness — Evan
- [ ] Freeze target provenance, development/locked partitions, metrics, nulls, and tuning budget before fitting — Susanna, supported by Evan
- [ ] Confirm PBTP EPIC generation, approved lineage manifest, genomic reference availability, and permitted output scope — Susanna and data custodian
- [ ] Assign figure/demo ownership during the first team check-in — Evan and Susanna

## Milestones

- **Day 1:** Verify data, freeze cohort/feature/protocol artifacts, pass software smoke tests, and produce the first development cancer-held-out prediction and null comparison.
- **Day 2:** Finish development LOCO, select and freeze the model, run locked CNS and eligible PBTP evaluation, and attempt CNV/fusion only if core results are stable.
- **Day 3:** Complete uncertainty, OOD, interpretation, figures, Shiny demo, documentation, presentation, and rehearsal; stop major model development at noon.

## Definition of Done

The project is complete when every reported result can be regenerated from versioned code and documented inputs; patient groups and outer cancers are isolated correctly; all learned preprocessing is fit within training data; the baseline is compared with a null; limitations and target provenance are visible; and the demo consumes approved precomputed outputs without protected identifiers.

If completed early, add the total-CN benchmark and fusion comparison on identical paired patients/splits. Only after those results are stable may the team run one predeclared train-only augmentation sensitivity using untouched real validation.

## Risks and Questions

- The full 41.5-GB historical matrix was not found inside the repository workspace during the 2026-09-13 audit; confirm the team's downloaded location and checksum before extraction.
- PBTP may use EPIC v2 rather than EPIC v1, which requires a separately versioned feature bridge.
- Published HRD labels may have incomplete upstream caller/build/parameter provenance and assay-specific error.
- Strong pooled accuracy may reflect cancer lineage, purity, batch, or platform rather than transferable HRD signal.
- R and required packages must pass smoke tests before biological fitting; the previous audit host had no working R runtime.
- Evan Savage (@esavage111) is the primary mentor/support contact; methodology and data-governance questions should be raised in the Team 12 Slack channel.

## Decision Rules

- Keep CpG elastic net unless another method improves development cancer-held-out macro-MAE consistently on identical splits without worse calibration/abstention.
- Do not inspect locked CNS or PBTP outcomes while selecting features, hyperparameters, architecture, thresholds, calibration, or OOD rules.
- Drop CNV/fusion if paired technical agreement or development incremental value fails the predeclared gate.
- Never call total-CN surrogates canonical HRDsum, HRD-LOH, or TAI.
- Never use synthetic samples as independent validation.
