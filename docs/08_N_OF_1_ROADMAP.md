# n-of-1 roadmap and output contract

The MVP consumes approved **processed** methylation features; production raw-EPIC inference requires a separately validated assay bridge. A new sample must not change transformations, features, normals, weights, calibration or OOD references.

Freeze a bundle containing model/version, training cohort hash, patient partitions, endpoint/source/build, exact probe list/order, technical masks, array manifests, normalization/reference files, imputation/scaling, coefficients, calibration residual rule and calibration membership, OOD reference, QC thresholds, software lock/session and code revision. Protect model artifacts until their release is reviewed.

Return: project sample ID; model version; predicted_reference_HRDsum; individual interval and validity note; QC status; fraction missing; OOD score/status; reason for abstention. Predicted LOH/LST/TAI appear only if separately trained and validated. HRD-high probability is unavailable until a separately calibrated reference-defined classifier exists. CNV plots/features require actual intensity-derived CN. Top coefficient contributions are additive model terms, not causal mechanisms.

Two uncertainty axes matter: estimator uncertainty and reference/assay error. Split conformal can quantify empirical residual spread under exchangeability; it cannot guarantee pediatric coverage. OOD should combine raw assay checks, frozen PCA distance, orthogonal reconstruction error, feature-range/missingness checks and subtype-aware retrospective error audit. Initial code's standardized-distance score is only an engineering heuristic. Set abstention thresholds using development/technical controls; quantify risk-versus-coverage on untouched data. Low genomic purity or unreliable CN signal can invalidate the respective branch while methylation QC may still pass; do not force fusion when one view fails.

Deployment stages: (1) processed public profiles; (2) paired raw 450K/EPIC reproducibility; (3) approved internal PBTP inference; (4) independently validated pediatric extension; (5) packaged prospective research tool. No claim of treatment selection without direct clinical/functional validation.

Shiny MVP in `app/app.R` reads precomputed results only, labels synthetic fixtures prominently, and shows prediction/interval/QC/OOD and actual-vs-predicted scatter where available. Do not allow unrestricted raw genomic upload for the presentation. Add CN chromosome plots only if approved real profiles are ready. Timebox UI to two hours; static figures are a complete fallback.
