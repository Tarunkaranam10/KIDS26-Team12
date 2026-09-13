# Methods decision matrix

All ratings are **recommendations/inferences**, not measured KIDS26 performance. Expected performance cannot be honestly ranked numerically before data. H/M/L indicate relative favorable suitability; leakage column describes the risk rather than a score.

| Candidate | Expected utility | Interpretability | Sample efficiency | Batch/450K-EPIC robustness | Cross-cancer expectation | Compute / implementation | n-of-1 | Leakage risk | Biology / three-day decision |
|---|---|---|---|---|---|---|---|---|---|
| A CpG-only | Strong baseline hypothesis | M-H, coefficients correlated | H with regularization | M after frozen manifest/bridge | Uncertain; tissue dominates | Low / easy after data | H if frozen | Selection/scaling outside folds | Primary |
| B methyl-CNV-only | May capture broad instability | H | H, few features | M; normals/purity dependent | Uncertain, more structural | Caller cost medium; regression low | H with fixed normals | Reference CN used as input; manual shifts | Secondary 1, gated |
| C fusion | Possible incremental information | M | M | Two sources of assay shift | No assumed improvement | Medium | H with both fixed views | Different subsets/late selection | Secondary 2, gated |
| D components | Reveals which scars predictable | H output interpretation | M; no extra independent N | Same as inputs | Uncertain | 3 fits | H | Error covariance/same-label overclaim | Post-core diagnostic |
| E RNA elastic net/ssGSEA | Orthogonal modality precedent | M | H for regularized model | Not 450K/EPIC; RNA processing mismatch | Adult evidence; pediatric unknown | Extra acquisition/integration | Conditional on fixed ranks/gene universe | Cohort normalization/source overlap | Defer unless ready |
| F PCA/AE/VAE | Compression; may encode tissue | L-M | PCA H; VAE lower | Not automatic harmonization | No guarantee of unseen-tissue transfer | PCA low; VAE high | Frozen encoder possible | Unsupervised full-cohort fitting | PCA sensitivity later; VAE post-event |
| G RF/XGBoost | Possible nonlinear gain on reduced data | M-L | M | No intrinsic protection | Extrapolation poor | Medium, extra tuning | H | Global feature selection/tuning | Not core |
| H elastic net | Stable p>>n baseline; extrapolates linearly | H subject to correlated CpGs | H | Must freeze transformations | Test, do not assume | Low; R glmnet | H | cv.glmnet alone does not nest preprocessing | Chosen estimator for A/B/C |
| I deep NN | Unproven benefit for this N/endpoint | L | L | Needs explicit domain work | High uncertainty | High | Possible | Large tuning/search surface | Cut from hackathon |

## Operational comparison

Beta scale is primary; M=log2((beta+epsilon)/(1-beta+epsilon)) is a predeclared sensitivity analysis only, with fixed epsilon and no endpoint-driven choice on heldout cancers. No global batch correction or outcome-based CpG filtering before splits. Technical masks can be fixed externally because they do not learn from these outcomes.

Train-only variance cap (default 5000, fixed before comparison) controls memory; compare a larger cap only in development tuning if time permits. The same selected-feature procedure runs independently inside each inner split. Elastic net alpha grid 0.1,0.5,1 and absolute lambda grid 0.01,0.1,1,10,100 are an initial engineering configuration, not a claim of optimality. Grid expansion belongs to development only. Log outcome transform is optional later; evaluate back-transformation bias.

Nulls: training mean, cancer-type means for within-known-type diagnostic comparisons (unseen type falls back to training mean), purity-only where observed, FGA-only for CN studies, and outcome permutations within cancer. Tumor type is not a primary predictor; unknown-type one-hot cannot represent novel cancers. Weighting cancers equally may improve transfer but trades pooled accuracy; choose that variant only inside development.

Decision rule: freeze the CpG baseline regardless of whether complex branches run. Compare B/C on identical eligible paired patients. A predeclared 5% relative macro-MAE improvement is a pragmatic promotion target, not a literature-established threshold; require consistency across cancers and report bootstrap uncertainty rather than claiming significance from tiny N.

Evidence basis: local L01/L42 constrain direct scar reconstruction; L37 supports association but targets MLPA; L10 supports assay/threshold caution; L43/L54 show single-sample methylation classifiers need not be deep; L64 does not validate HRD augmentation. See the evidence CSV for gaps.
