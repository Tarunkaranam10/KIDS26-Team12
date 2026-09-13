# Risk register

| Risk | Detection | Mitigation | Fallback |
|---|---|---|---|
| Assay/aliquot mismatch | Duplicate/specimen/timepoint audit | Explicit hierarchy, unique join, retain IDs | Exclude unresolved |
| Missing HRD labels | Cohort flow and label sum checks | Published target plus source audit | Smaller labeled cohort |
| IDAT access | Live metadata, both channels, access status | Selected public subset + fixed normals | CpG-only |
| Poor methylation CN | Paired length-weighted/arm/breakpoint agreement | Blind caller pilot, purity strata | Drop CN branch |
| Allelic information inadequate | Copy-neutral LOH missed, sparse informative loci | Matched-normal SNP experiment | Predict LOH/TAI; no direct calculation |
| 450K/EPIC shift | Paired processing disagreement | Exact manifest intersection, frozen per-sample processing | Retrospective 450K claim only |
| Adult/pediatric shift | Held-out CNS/PBTP calibration and OOD | Locked external testing, abstention | Negative transfer result |
| Tumor-type confounding | Strong pooled but weak within-cancer signal | LOCO macro metrics, type/null controls | Narrow disease-specific hypothesis |
| Purity/ploidy | Error vs independent purity/ploidy | QC/dilution and stratified reporting | Abstain or restrict use |
| Batch effects | Errors by center/processing/platform | Fixed references, technical replicates | Separate data tracks |
| Leakage | Group overlap, selected test features, generator donors | Nested preprocessing and group tests | Rerun corrected pipeline; invalidate prior metrics |
| Overfitting | Outer-negative R2, unstable features | Small fixed search, regularization | Mean/baseline comparison only |
| Insufficient PBTP N | Independent donor count | Cluster bootstrap, descriptive results | Defer generalization claim |
| Synthetic artifacts | Bad covariance, near-duplicates, no held-out gain | Real-only primary; train-fold generation | No augmentation |
| Miscalibration | Coverage/slope/PR prevalence audit | Separate calibration, domain warning | No probability/interval claims |
| UI consumes event | Missed scientific milestones | Precomputed outputs, two-hour budget | Static plots |
| Scar truth is historical | Discordance with function/response | Orthogonal assays and wording | Restrict to scar-burden prediction |
| Unavailable R/runtime | R tests/session capture fail | Pre-event validated R environment | Acquisition/engineering only |
| Public data exposure | Git status/staged file review | Ignored data/results, external protected storage | Approved synthetic demo only |
