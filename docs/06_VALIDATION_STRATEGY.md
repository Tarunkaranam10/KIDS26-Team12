# Validation, harmonization and analysis contract

## Partitions and estimand

Primary estimand: prediction error for continuous, reference-defined unadjusted HRDsum in a new cancer type. Secondary: within-cancer error and association; reference HRD>=42 classification is exploratory, not pediatric clinical HRD. Training is non-CNS TCGA. Development validation is nested LOCO within non-CNS cancers. GBM/LGG are locked domain tests. PBTP is a final locked external test, architecture and preprocessing frozen first. Any adaptation using PBTP features, labels or test-driven QC thresholds creates an exploratory/development cohort and requires new independent testing.

Before fitting, save immutable patient/specimen groups and assay eligibility. Group all aliquots of a donor. For PBTP group all derived lineages; primary metric uses patient tumors; model lines are secondary and clustered. No synthetic patients or masks derived from held-out donors in training. Neither supervised nor unsupervised transforms may learn from held-out data under the strict transfer claim.

Use outer leave-one-cancer-out. Tune on cancer-grouped inner folds where enough cancers exist; otherwise grouped-patient folds with explicit reduced transfer evidence. In every inner training split fit missingness/variance selection, medians, scales, latent encoders and any augmentation anew. `cv.glmnet` after global feature selection is insufficient. Select hyperparameters on mean per-cancer MAE; report all outer cancers, including failures. Do not repeatedly use outer holdout results to tune a large architecture search and call them unbiased final performance.

Training-only mean is the primary null; show MAE/RMSE improvement versus null per cancer. Known-type means, purity-only, FGA-only and within-cancer label permutations diagnose confounding. Pooled Pearson correlation can be dominated by between-tissue differences. Report macro-averaged cancer metrics and within-cancer correlations; macro R2 excludes undefined strata but reports their count. Do not discard low-variance/negative-R2 cancers as outliers.

## Metrics and uncertainty

For every run: cohort/sample flow, training/test N by cancer, independent patient N, input/selected feature count, target distribution, preprocessing provenance, hyperparameters and fold assignments. Continuous: MAE, RMSE, R2 (undefined when outcome variance zero), Pearson/Spearman with valid-pair counts, calibration intercept/slope, residual bias, interval coverage/width and abstention fraction. Compare methods paired by patient. Bootstrap patients within cancer; PBTP bootstrap donor lineages. Tiny-N CIs are descriptive and may be unstable; do not invent a power calculation without anticipated effect and variance.

Binary secondary: ROC-AUC, PR-AUC with prevalence baseline, sensitivity/specificity, PPV/NPV, Brier score, reliability plot and threshold sensitivity. Metrics requiring both classes are NA in one-class strata. A regression prediction is not a calibrated HRD-high probability. Fit any probability mapping solely on separate training calibration data and validate it; omit it from MVP outputs if unavailable. Never infer probabilities from a point estimate or an interval alone.

Use split conformal residual intervals on a reserved calibration group, untouched by fitting and tuning. Keep that set excluded from the final fitted model or recalibrate on new independent data after refitting. The finite-sample quantile uses ceil((ncal+1)*(1-alpha)); insufficient calibration N produces an unbounded/unavailable interval, not an unjustifiably finite one. Exchangeability assumptions do not hold automatically for unseen tissue or adult-to-pediatric transfer; coverage there is measured, not guaranteed. A confidence interval for mean prediction is not an individual prediction interval.

## Methylation preprocessing / 450K to EPIC

Confirm EPIC generation explicitly; v2 has duplicated/revised probe IDs and mappings. Freeze common probe IDs using authoritative manifests, coordinates, strand and build. Avoid using an imputed or stripped ID to merge distinct probes. Primary features are autosomal CpGs; remove technical SNP, cross-reactive/non-unique and manifest-masked probes using versioned technical annotations. Use ancestry-aware common-variant annotations rather than treating one population as universal. Preserve genotype-sensitive channels in a separate experimental branch.

Raw production candidate: one pinned SeSAMe per-sample pipeline with documented background, dye-bias, detection and masking steps; fixed conumee2 normals for CNV. Benchmark minfi single-sample noob as a technical comparator if needed. Freeze all reference files and package versions. Do not use target-cohort ComBat, joint quantile normalization, refitted PCA or new-cohort centering at deployment. Batch correction learned on historical cohort betas is not automatically exportable to a raw EPIC sample.

Beta values are primary for bounded scale and interpretability. Train-fold M-value sensitivity uses fixed epsilon only. Fit missingness filter and imputation medians on training; maintain feature order and reject a new sample with too many absent model features (initial engineering limit 5%, to be calibrated on development technical data). Missingness can itself encode deletions/platform; report it, do not let imputation hide assay failure. Technical sample failure thresholds must be predeclared from caller/reference guidance and a development QC pilot; published quality flags alone are incomplete.

Purity/ploidy are initially QC/confounder variables, not required model predictors because the deployed sample may lack independent genomic purity. Never adjust the target to remove a biological signal without defining a new estimand. Investigate errors across purity/ploidy strata and using independently measured purity where possible. Masking high-CNV regions can be a development sensitivity analysis to distinguish epigenetic from CN-mediated association.

## Frozen inference

Input manifest -> per-sample QC -> fixed transformation -> exact feature mapping -> stored training medians/scales -> saved model -> calibration/OOD -> structured report. Permuting input feature order or scoring one sample alone vs in a batch must give the same result. Reject duplicates, incompatible build/platform, impossible beta values and excessive missingness. A processed-beta MVP is not yet a validated raw-IDAT calculator.

The initial R code provides nested elastic-net fitting and a conservative standardized-distance OOD heuristic. It does not implement raw EPIC preprocessing, calibrated binary probabilities, or a validated domain-shift detector. These are explicitly gated extensions. OOD thresholds themselves must be assessed on held-out technical/domain samples; passing OOD does not prove reliability.
