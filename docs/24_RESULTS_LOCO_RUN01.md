# 24 — Results: LOCO run 01 (2026-09-17)

**Run.** LSF array `323078995`, 30/30 folds DONE, merged 06:20 CDT by
`scripts/watch_loco_array.sh`. Development cohort n = 7,065 across 30 non-CNS
cancer types. Locked CNS partition (GBM+LGG, 642 samples) **still unopened.**

**Verdict.** The tissue confound is **real but not disqualifying.** Modelling
HRDsum continues, with a changed claim: this is a **within-tissue relative
ranker**, not an absolute HRD calculator. See §6 for what must change.

---

## 1. Headline numbers

| Quantity | Value | File |
|---|---|---|
| Pooled MAE (model) | **9.049** | `pooled_null_panel.tsv` |
| Pooled MAE (tissue-mean null) | 9.890 | `pooled_null_panel.tsv` |
| **Skill vs tissue-mean null** | **+0.085** | `pooled_null_panel.tsv` |
| Skill vs training-mean null | +0.328 | `pooled_null_panel.tsv` |
| Macro MAE (unweighted over tissues) | 8.566 | `macro_metrics.txt` |
| **Within-tissue Pearson** | **0.612** | `pooled_within_tissue_metrics.tsv` |
| Within-tissue Spearman | 0.601 | `pooled_within_tissue_metrics.tsv` |
| Within-tissue MAE | 7.881 | `pooled_within_tissue_metrics.tsv` |
| Within-tissue permutation p | **0.001** (1000 perms) | `pooled_within_tissue_permutation.tsv` |
| Lambda at path boundary | **0 of 30 folds** | `loco_metrics.tsv` |
| Reportable rate | 99.2% | `loco_predictions.tsv` |

---

## 2. The apparent contradiction, and its resolution

Two numbers look incompatible:

- Skill vs tissue-mean null = **+0.085** — barely better than guessing each
  cancer type's average.
- Within-tissue correlation = **0.612**, permutation p = 0.001 — strongly
  better than chance at ranking patients inside a tissue.

Both are correct. They disagree because the model gets the **ranking** right and
the **absolute level** wrong.

Per-tissue calibration offset (`bias = mean(predicted − actual)`):

| Statistic | Value |
|---|---|
| Mean absolute offset | **3.46 HRD units** |
| SD of offset across tissues | 4.83 |
| Range | −15.44 (KIRC-like low) to +8.29 (PCPG) |
| Prediction shrinkage (mean sd_pred/sd_actual) | 0.772 |

Remove that per-tissue offset and pooled MAE falls **9.049 → 7.881**. The offset
alone costs ~1.17 MAE units — more than the entire margin over the tissue null.
So the model *has* learned HRD-relevant signal; it is squandered by systematic
mis-levelling of each tissue.

`results/loco_run01/figures/within_vs_absolute.png` shows this directly.

---

## 3. How much of the model is just tissue identity?

| Test | Result | Reading |
|---|---|---|
| Variance of **predictions** explained by tissue identity alone | R² = **0.562** | 56% of what the model says is recoverable from the label "this is a breast tumour" |
| Variance **not** explained by tissue | **0.438** | 44% is something else |
| Variance of **true HRDsum** explained by tissue identity | R² = **0.341** | The confound is real in the biology, not invented by the model |
| Within-tissue permutation test | p = 0.001 | Residual signal survives after tissue is removed |

The model leans on tissue harder than the truth warrants (0.562 vs 0.341) — it
**over-weights lineage**. But 44% of its behaviour is not tissue, and that
portion is what produces r = 0.61 within tissue. A pure tissue-lookup model
would show within-tissue r ≈ 0 and permutation p ≈ 0.5. It does not.

---

## 4. Per-tissue breakdown

29 of 30 tissues have within-tissue r > 0; 26 of 30 exceed r > 0.3; median
r = **0.551**. Full table: `per_tissue_within_correlation.tsv`.

**Best:** KICH 0.803, UCEC 0.758, STAD 0.738, LUAD 0.721, BLCA 0.675, PRAD 0.663.

**Failures:** THCA **−0.032** (bias +8.05), PCPG 0.110 (bias +8.29),
CHOL 0.161 (n=35), TGCT 0.162.

The failure pattern is coherent: THCA and PCPG are genomically quiet, low-HRD
tumours, and the model **over-predicts both by ~8 units**. Trained mostly on
scarred tumours, it cannot express "this genome is calm." Combined with the
zero-floor problem (§5) this is the single most actionable defect.

---

## 5. Limitations that constrain interpretation

**Ovarian is effectively missing — n = 10.** OV is the canonical HRD cancer and
the primary clinical use case for HRD testing. TCGA ovarian methylation is
predominantly 27k-array, so it was excluded by the 450k probe bridge. Any claim
about HRD-high tumours rests on UCEC/BRCA/STAD instead. **This must be stated in
the presentation.** Similarly thin: CHOL 35, DLBC 47, UCS 56.

**HRDsum is bounded at zero and 14.3% of samples sit exactly at 0.** The elastic
net is unbounded and emits negative predictions (visible in the figure). MAE is
inflated by predictions that are impossible a priori. Clipping at 0, or modelling
`log1p(HRDsum)`, is a free improvement not yet applied.

**Purity is a live, unresolved concern.**

| Quantity | Value |
|---|---|
| cor(prediction, purity) within tissue | **0.165** |
| cor(true HRDsum, purity) within tissue | **0.019** |

Predictions track tumour purity roughly **8× more strongly than the truth does**.
Some of what the model reads is "how much tumour is in this sample," not HRD.
Skill degrades monotonically as purity rises:

| Purity stratum | n | MAE model | MAE tissue null | Skill |
|---|---|---|---|---|
| Low (0.08–0.51) | 2351 | 8.599 | 10.611 | **+0.190** |
| Mid (0.52–0.73) | 2344 | 9.184 | 9.726 | +0.056 |
| High (0.74–1.00) | 2202 | 9.484 | 9.112 | **−0.041** |

Note the direction: model MAE *worsens* with purity (8.60→9.48) while the tissue
null *improves* (10.61→9.11). In the cleanest, highest-purity samples the model
**loses to the tissue mean.** That is the opposite of what a genuine biological
signal should do and is not yet explained. Purity-matched subset
(0.5–0.8, n=3312) retains skill +0.045, so the effect is attenuated but not
abolished.

**No CNS validation yet.** The locked partition is untouched. Everything above
is development-cohort performance and does not establish transfer.

---

## 6. Decision and required changes

**Continue modelling HRDsum.** Justification: within-tissue r = 0.61 with
permutation p = 0.001 across 29/30 tissues is not a tissue-lookup artefact, and
0/30 folds hit a lambda boundary so the fit is not grid-limited.

Required before any result is presented or the CNS lock is opened:

1. **Reframe the claim.** Report as a within-tissue relative ranker. Lead with
   within-tissue r = 0.61 and permutation p = 0.001, not with skill = +0.085.
   Do **not** describe the output as an HRD score comparable across tissues.
2. **Fix per-tissue calibration.** Worth ~1.17 MAE units. Must be fitted inside
   the LOCO loop on training tissues only — a per-tissue offset fitted on the
   held-out tissue is leakage and invalidates the fold.
3. **Enforce the zero floor.** Clip at 0 or model `log1p(HRDsum)`.
4. **Resolve the purity inversion** (§5). Until explained, the high-purity
   result blocks any strong biological claim.
5. **State the OV limitation** in every presentation of these numbers.
6. **Keep the CNS lock closed** until 1–4 are done. It opens once.

---

## 7. Provenance

All numbers from `results/loco_run01/`, array `323078995`, commit `696e05e`
plus this run's outputs. Matrix `data/processed/beta.tsv`
(sha256 `e3642d30…`, `engineering_only: false`), 336,480 probes × 7,707 samples.
Per-fold artefacts in `results/loco_run01/folds/` (30 × 5 files).
