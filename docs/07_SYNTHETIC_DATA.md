# Synthetic data and mechanistic simulation

**Decision:** no synthetic augmentation in the primary model. First complete a real-only baseline and unseen-cancer evaluation. The simplest optional augmentation is constrained interpolation in a training-fitted PCA space, not VAE pretraining.

| Method | Relevant evidence / limitation | Decision |
|---|---|---|
| SMOTE | Local L62 binary tumor/normal classification; endpoint/nesting do not establish HRD regression utility | No raw-CpG SMOTE primary; compare weighting/resampling first |
| Bootstrap | Quantifies training instability; resampled donors are not independent new patients | Use for uncertainty, not increased N claims |
| Technical noise/masking | Local Sturgeon/TUCAN simulate assay sparsity, with real validation | Useful robustness tests with source-donor partitioning |
| PCA latent interpolation | Simple geometric regularization; assumes locally smooth outcome | Optional post-core sensitivity |
| AE/VAE | Representation learning feasible, may encode tissue more than HRD | Post-event |
| Conditional VAE | Conditions on observed training labels; requires enough data per condition | Post-event benchmark |
| GAN/conditional GAN | OncoGAN generates factorized genomic features, not methylation-HRD patients | No three-day core dependency |
| Diffusion | DiffuCpG is inpainting; good reconstruction does not prove HRD utility | Defer |
| Methylation-specific generators | Must preserve CpG covariance, genomic context and target relation | Benchmark on real heldout data before adoption |

[ SyntheVAEiser ](https://link.springer.com/article/10.1186/s13059-024-03431-3) concerns VAE-generated **gene expression** for cancer-subtype prediction, not validated HRD methylation augmentation. Its leave-cancer experiment also fine-tunes on 40 real target-cancer cases; this is not zero-shot transfer to untouched PBTP. This weakens the previous plan's direct analogy. The Precious2GPT page and the 2026 Epi-HRD score page were not retrievable during the source check; do not treat their detailed methods as verified.

## Optional PCA experiment

Inside each inner training split only: technical allowlist -> imputation/scaling -> PCA with predeclared low rank (e.g. <=20 and <Ntrain) -> choose near neighbors within training cancer and near HRD values -> interpolate latent vectors and observed labels with the same coefficient. Local linearity of HRD along the interpolation is an assumption; do not relabel with the predictor under evaluation. Limit augmentation to <=0.5x real N as an initial experimental cap, not an established optimum. Compare no augmentation and sample weighting on identical folds. Do not generate PBTP-like patients using held-out PBTP profiles.

Fidelity: beta range after decoding, marginal distributions, CpG correlation and regional/chromosomal covariance, distances to donors, PCA/UMAP fit on training only, nearest-neighbor duplicate/memorization audit, HRD relationships, and subgroup representation. Measure train-synthetic/test-real only on untouched donors, alongside real-only vs real+synthetic on the same held-out cancers. Attractive embeddings and large synthetic N are not biological validation. Nearest-neighbor privacy checks do not prove privacy; internal synthetic data require governance review before release.

## Mechanistic simulator: a separate engineering product

`scripts/cnv_features.py` computes explicitly named total-CN summaries; it is not scarHRD. `scripts/simulate_scars.py` produces analytic segment fixtures and noisy methylation-like total-ratio observations. Include balanced diploid, copy-neutral LOH, telomeric allele-imbalanced regions, balanced gains, whole-chromosome LOH, centromere-spanning AI, and >=10 Mb state transitions. Ground-truth event labels in fixtures are hand-defined under stated simplified geometry and must not be used as a substitute for canonical package validation.

For mixture purity p, total depth ratio relative to diploid normal is (p*(A+B)+2*(1-p))/2. At a heterozygous locus, minor-allele fraction is (p*B+(1-p))/(p*(A+B)+2*(1-p)). Balanced 1/1 and copy-neutral 2/0 have identical total signal at every purity but different BAF. This proves non-identifiability from total intensity alone; adding more total-CN probes cannot recover the lost distinction.

Stress dimensions: purity 0.2/0.5/0.8/1, tumor ploidy/background, Gaussian signal noise, short false segments, missing probe intervals, breakpoint shifts and coarser bins. Quantify component error against pinned scarHRD on known allele-specific input, and against total-CN proxies after erasing alleles. Current simulator covers the identifiability fixtures and purity/noise/probe-grid degradation; full ploidy-aware canonical scar benchmarking remains post-core work. Do not equate an engineering test pass with clinical validity.

## Differential evidence update — 2026-09-10

New L59 reviews current synthetic cancer data and separates fidelity, downstream utility, and privacy; it identifies validation standardization as unresolved. L61 omicsGAN reports an internal methylation AUC increase in 173 TCGA AML cases, but has no HRD target, external cohort, pediatric transfer, or 450K/EPIC bridge. L65 demonstrates the correct comparison structure—real-only versus real-plus-synthetic training on an untouched real test set—but in histology rather than methylation. L63 adds general opportunities and risks without a directly transferable benchmark.

These papers do not justify SMOTE, VAE, conditional VAE, GAN, diffusion, or a methylation generator in the three-day critical path. The optional gated experiment remains train-only PCA-space interpolation, capped and fit within each training fold, compared with real-only training on exactly the same untouched real validation cancers. It may begin only after the real baseline, held-out-cancer result, null comparator, and transfer audit are complete. Synthetic cases are never an independent validation cohort.
