# R pipeline and Shiny MVP audit — 2026-09-10

## R implementation

Entry points are `scripts/train_baseline.R`, `scripts/predict_frozen.R`, and `tests/smoke_model.R`; fold-safe helpers are in `R/model.R`. Required R packages are `glmnet`, `data.table`, and `jsonlite`; `shiny` is required only for the presentation layer. `scripts/setup.R` initializes `renv`, installs these packages, and writes session information. A lockfile must be snapshotted only after installation and tests succeed on the intended analysis host.

| Concern | Code audit | Status |
|---|---|---|
| Outer validation | Each non-CNS cancer is excluded before `fit_en`; GBM/LGG are excluded from all development fits and are not evaluated by the training script | Correct by inspection |
| Inner tuning | With at least three training cancers, each inner fold is an entire cancer; alpha/lambda selection minimizes macro-MAE across held-out cancers | Correct by inspection |
| Missingness and imputation | Feature missingness and medians are learned inside each inner-training split, then refit on the complete outer-training set | No outer-test fitting |
| Variance filtering | Variance ranking/top-5,000 selection occurs in `fit_preprocess` on the current training split | No outer-test fitting |
| Scaling | Centers and scales are learned in `fit_preprocess`; `glmnet` receives `standardize=FALSE` | No double/global scaling |
| Null comparator | Training-set mean is saved in each outer bundle and applied to its held-out cancer | Leakage-safe |
| Final calibration | A fixed-seed, cancer-stratified 20% non-CNS reservation is excluded from final model fitting; its residuals set the conformal radius | Correct; coverage does not transfer automatically to PBTP |
| OOD | Training standardized-distance 99th percentile is stored and applied without refitting | Heuristic, explicitly not a calibrated domain-shift test |
| Frozen inference | Features are matched by name/order; missing features use stored medians; excess missingness/OOD withholds display; model MD5 and provenance are emitted | n-of-1 compatible at model layer |
| Coefficients | Final intercept and feature coefficients are exported to `frozen_coefficients.tsv` | Added in this audit |
| Reproducibility | Seed 260910, fold table, final partitions, tuning grid, session information, probe provenance, and serialized bundles are saved | Implemented; package lock pending runtime |

`tests/smoke_model.R` creates a 90-sample/30-CpG engineering matrix and checks preprocessing immutability, patient uniqueness, single-sample/batch prediction equality, feature-order invariance, conformal behavior, and metric edge cases. These tests validate software mechanics only. They do not validate biological prediction, 450K/EPIC transport, interval coverage under domain shift, or pediatric use.

## Execution result

Python compilation and 17 Python tests pass. A native R executable was not found in standard Windows, user-program, project, or conda locations. The installed WSL distribution also lacks R, and package installation requires interactive administrator authentication. Therefore R source and leakage boundaries were inspected, but `parse()`, package loading, the R smoke test, model training, and Shiny rendering could not be executed on this host. No biological model result is claimed.

The first R-host commands remain:

```sh
Rscript scripts/setup.R
Rscript tests/smoke_model.R
Rscript -e "renv::snapshot(prompt=FALSE)"
```

Do not start a biological fit unless all three succeed and the full matrix/probe provenance gate has passed.

## Shiny MVP

`app/app.R` reads only a precomputed TSV named by `KIDS26_DEMO_RESULTS`; otherwise it displays clearly labeled synthetic fixtures. It does not train, tune, upload, or write analytical data. The current MVP displays:

- selected sample and provenance;
- predicted reference HRDsum and interval when reportable;
- independent reference HRDsum and residual when supplied;
- QC/OOD warning and withheld prediction state;
- aggregate validation N/MAE/RMSE and per-cancer N/MAE when labels/cancer types are supplied;
- prediction-versus-reference scatter and predicted-score distribution.

The app omits CNV plots unless a precomputed, technically validated CN branch is later approved. Protected PBTP identifiers and raw data are not required for demonstration; public, synthetic, or explicitly approved de-identified aggregate outputs are sufficient. Runtime rendering remains unverified because R/Shiny are unavailable on this host.
