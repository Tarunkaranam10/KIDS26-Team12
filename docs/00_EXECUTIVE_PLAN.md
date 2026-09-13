# KIDS26: executable scientific plan

## Decision

Build a **research predictor of reference genomic-scar burden**: shared-probe methylation -> elastic-net continuous HRDsum regression. Keep methylation-derived CNV-only and methylation+CNV elastic-net fusion as the only two secondary model families for the hackathon, conditional on paired intensity/CNV data. Do not call this a direct HRD calculator, current functional HR assay, or treatment-response test.

**Documented facts:** local Abkevich defines HRD-LOH using allele-aware LOH regions; CopyNumber450kCancer explicitly notes absent BAF; the closest local methylation study (Chen) targets MLPA BRCAness probability rather than canonical HRDsum. Rempel's full TCGA scar analysis uses SNP-array ASCAT, with WES in an OV comparison subset. See [evidence matrix](literature_evidence.csv) and [paper digest](01_LITERATURE_EVIDENCE.md).

**Inference:** methylation may predict scars, but tissue identity, purity, proliferation and platform can explain an apparently strong pooled association. A good random-split correlation is insufficient.

**Recommendation:** acquire independent labels first, test unseen cancer transfer, retain a no-CNS development boundary, and freeze before PBTP. A failed transfer experiment is a useful scientific result if reproducible and adequately powered.

## Primary and secondary experiments

1. CpG beta -> continuous HRDsum, elastic net, no cancer-type predictor. Start with technically allowed shared autosomal probes, training-fold variance cap, imputation and scaling. Do not require raw intensities for this retrospective MVP.
2. conumee2 total-CN structural features -> HRDsum, elastic net. Go only after blinded reference agreement/QC benchmark. SNP6/ABSOLUTE features are reference-only; feeding them as predictors would defeat the methylation-only input objective.
3. Fusion on exactly the same paired samples and splits. Promote only if development LOCO macro-MAE improves consistently without worse failure/calibration behavior; otherwise preserve CpG-only. Final locked tests do not select the winner.

Component regressions, RNA, nonlinear trees, PCA augmentation, VAE and diffusion are lower priority diagnostics/extensions, not equal-priority deliverables. Uncertainty and OOD are part of the reporting contract, but unsupported clinical probabilities must be unavailable rather than fabricated.

## Architecture and status

Python standard library handles acquisition, provenance, joins and engineering tests. R/Bioconductor handles methylation/CNV and glmnet modeling; Shiny consumes precomputed approved outputs. Avoid a second ML stack. The initial implementation includes data acquisition and preparation, a nested R baseline, tests, a total-CN feature utility and an engineering simulator, and a lightweight demo.

There is no fitted biological model or PBTP result yet. Differential reconciliation now covers 60 local PDFs: 35 from the prior review, 24 newly available entries from the original 69-row map, and one additional chordoma paper. The evidence matrix has 70 rows and preserves unavailable methods and supplemental gaps. See [the differential reconciliation](15_NEW_LITERATURE_RECONCILIATION.md). The historical .doc/.docx takeaways are planning hypotheses; the original CSV remains the bibliographic source for L01-L69.

## 2026-09-10 evidence delta

The primary architecture did **not** change. New allele-aware TAI, HRDsum, intratumor, and GSS methods strengthen the conclusion that canonical LOH and TAI require allelic information. New laboratory consensus guidance and CNS classifier methods **partially refine** implementation: target caller/build/parameters and tumor-specific thresholds are mandatory provenance; the 450K/EPIC candidate mask is versioned; feature coverage and abstention are explicit; and patient/model derivatives remain grouped. No new synthetic-omics study supplies sufficient external HRD transfer evidence to move augmentation into the three-day critical path.

## Meaningful success

- **Minimum:** reproducible matched public dataset, auditable labels, leakage-safe baseline and at least one real held-out cancer evaluation with null controls; report failure if it fails. A six-sample smoke test alone does not meet this threshold.
- **Strong:** multi-cancer transfer audit, frozen eligible PBTP patient-tumor evaluation, interval/abstention diagnostics, and clear limitations.
- **Stretch:** validated methylation-CNV incremental benefit, component decomposition, or a predeclared augmentation improvement on untouched real cases.

## Prior plan critique

| Prior choice | Assessment and change | Main risk | Falsification / simpler fallback |
|---|---|---|---|
| Direct CpG -> continuous sum | Keep as hypothesis; direct local precedent has a different MLPA target | Tissue and purity confounding | Fails within-cancer/LOCO versus training-only mean; publish negative transfer result |
| CNV-derived canonical scars | Reject canonical LOH/TAI from total CN; retain total-CN surrogate features | Non-identifiability, centromeres and resolution | Paired CN benchmark fails; drop CN branch |
| Fusion as likely final model | Remove presumption of superiority | Paired-subset selection and feature scale | No stable development LOCO gain; keep CpG-only |
| Predict components then sum | Diagnostic extension only | Correlated errors; same labels do not add data | Worse sum MAE or nonsensical components; retain direct sum |
| RNA comparator | Useful orthogonal modality; not independent truth | Same patients/labels; source overlap; cohort normalization | No compatible processed RNA; defer |
| VAE first synthetic approach | Reject for three-day core | Insufficient N, leakage, tissue encoding | No real-only baseline benefit; no augmentation |
| LOCO | Keep, strengthen with nested cancer-aware tuning and macro metrics | Hyperparameter reuse across outer holdouts | Pooled metric only; require per-cancer results |
| GBM/LGG holdout then PBTP | Keep locked hierarchy | Looking at CNS and revising while still calling it locked | If adapted, relabel CNS development and retain untouched PBTP |
| Frozen single-sample pipeline | Keep; historical betas do not prove deployability | Normalization/platform mismatch | No bridge validation; demo processed public profiles only |
| 42 probability and intervals | Secondary only, explicitly reference-defined | Threshold portability and false uncertainty | No calibration data; return unavailable |
| PBTP labels available | Prior handoff asserted availability; user says potentially available | Missing matched tumors/genomics | Internal pre-event manifest gate; no invented PBTP N |
| Generator sees heldout features unsupervised | Disallow for strict transfer claim | Transductive leakage even without labels | Fit every learned step on training only |

The root Word map's suggestion that methylation CNV leads directly to HRDsum omits the allelic requirement. The numbered markdown plan largely recognizes it already. The main improvements are sharper target provenance, narrower scope, strict preprocessing boundaries, and falsifiable transfer tests.
