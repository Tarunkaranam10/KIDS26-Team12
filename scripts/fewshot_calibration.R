# =============================================================================
# scripts/fewshot_calibration.R - Does a handful of labelled samples from an
# unseen tissue recover its calibration offset? (C1 follow-up)
# =============================================================================
#
#   Rscript scripts/fewshot_calibration.R        (from the repository root)
#
# BACKGROUND
#   docs/25 established that the per-tissue offset (mean 3.46 HRD units) cannot
#   be predicted from LABEL-FREE tissue covariates: every configuration lost to
#   a constant, best leave-one-tissue-out R^2 = -0.116.
#
#   This script tests the remaining practical option. If we are willing to
#   measure HRD on k samples from a new tumour type, does that buy us enough to
#   correct the offset for everyone else of that type?
#
# WHY THIS NEEDS NO REFIT
#   A per-tissue offset correction is a post-hoc shift of predictions that
#   already exist. We can therefore answer this exactly from run 01's held-out
#   predictions, with no cluster time at all. Every prediction used here was
#   made by a model that never saw the tissue it is predicting, so the LOCO
#   guarantee is intact.
#
# -----------------------------------------------------------------------------
# THE COMPARATOR PROBLEM - THE CRUX OF THIS ANALYSIS
# -----------------------------------------------------------------------------
# It is tempting to compare few-shot-corrected predictions against the ORIGINAL
# tissue-mean null and declare victory. That comparison is rigged, because it
# gives the model k labels while giving the null none.
#
# If a clinician has k labelled samples from a new tumour type, they can ignore
# our model entirely and just predict the mean of those k labels for every
# future patient. THAT is the honest competitor, and it gets better with k at
# the same time our correction does.
#
# So this script reports both, and the decision rests on the fair one:
#
#   model_uncorrected   run 01 behaviour, no labels used
#   model_fewshot       offset estimated from k labels, applied to the rest
#   model_oracle        offset estimated from ALL of the tissue's labels
#                       (an upper bound nobody can achieve in practice)
#   null_fewshot        predict mean(k labels) for everyone   <- FAIR COMPETITOR
#   null_oracle         predict the tissue's true mean for everyone
#
# Everything is evaluated ONLY on the held-out remainder of the tissue, never on
# the k samples used to estimate the offset, because those are now "training".
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

pred_path <- "results/loco_run01/loco_predictions.tsv"
if (!file.exists(pred_path)) stop("Run 01 predictions not found: ", pred_path)
p <- read.delim(pred_path)

set.seed(8123)
K_VALUES <- c(3L, 5L, 10L, 20L)
B        <- 200L   # random draws per (tissue, k); the spread across draws is
                   # itself a headline result, not just noise to average away
MIN_EVAL <- 10L    # need enough held-out samples for a stable MAE

# -----------------------------------------------------------------------------
# One (tissue, k, draw): pick k patients, estimate the offset from them alone,
# correct everyone else, and score on those others.
# -----------------------------------------------------------------------------
one_draw <- function(df, k) {
  idx  <- sample.int(nrow(df), k)
  shot <- df[idx, , drop = FALSE]     # the k "labelled" samples
  ev   <- df[-idx, , drop = FALSE]    # everyone else - the evaluation set

  # The offset a practitioner could actually compute from k labelled samples.
  off_hat <- mean(shot$predicted_reference_HRDsum - shot$actual)
  # The offset they could compute with the whole tissue labelled (unattainable).
  off_oracle <- mean(df$predicted_reference_HRDsum - df$actual)

  c(
    model_uncorrected = mean(abs(ev$actual - ev$predicted_reference_HRDsum)),
    model_fewshot     = mean(abs(ev$actual - (ev$predicted_reference_HRDsum - off_hat))),
    model_oracle      = mean(abs(ev$actual - (ev$predicted_reference_HRDsum - off_oracle))),
    # The fair competitor: k labels, no model.
    null_fewshot      = mean(abs(ev$actual - mean(shot$actual))),
    null_oracle       = mean(abs(ev$actual - mean(df$actual))),
    # How far the k-sample offset estimate lands from the true tissue offset.
    offset_error      = abs(off_hat - off_oracle)
  )
}

res <- list()
for (t in sort(unique(p$cancer_type))) {
  df <- p[p$cancer_type == t, , drop = FALSE]
  for (k in K_VALUES) {
    if (nrow(df) - k < MIN_EVAL) next   # tissue too small to test this k
    draws <- replicate(B, one_draw(df, k))
    res[[length(res) + 1]] <- data.frame(
      cancer_type = t, k = k, n = nrow(df),
      t(rowMeans(draws)),
      # Spread of the offset estimate across draws: the instability a
      # practitioner would actually experience with an unlucky sample.
      offset_error_sd = sd(draws["offset_error", ]),
      stringsAsFactors = FALSE
    )
  }
}
res <- do.call(rbind, res)

dir.create("results/fewshot", showWarnings = FALSE, recursive = TRUE)
write.table(res, "results/fewshot/fewshot_by_tissue.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# -----------------------------------------------------------------------------
# Aggregate. Weight tissues equally rather than by n: the question is "does this
# work for a NEW tissue", and each tissue is one trial of that question. Sample
# weighting would let BRCA and THCA dominate the answer.
# -----------------------------------------------------------------------------
cat("\n=== FEW-SHOT CALIBRATION (mean over tissues, equal weight) ===\n")
# Per-tissue win flags must be computed BEFORE aggregating. Computing them after
# summarise() would collapse 29 comparisons into one and always report 0 or 1.
res$beats_fair  <- res$model_fewshot < res$null_fewshot
res$helps_model <- res$model_fewshot < res$model_uncorrected

agg <- res |>
  group_by(k) |>
  summarise(
    tissues            = n(),
    model_uncorrected  = mean(model_uncorrected),
    model_fewshot      = mean(model_fewshot),
    model_oracle       = mean(model_oracle),
    null_fewshot       = mean(null_fewshot),
    offset_err         = mean(offset_error),
    offset_err_sd      = mean(offset_error_sd),
    # How many individual tissues benefit - the mean can hide a split decision.
    n_helps_model      = sum(helps_model),
    n_beats_fair_null  = sum(beats_fair),
    .groups = "drop"
  ) |> as.data.frame()
print(agg, digits = 3, row.names = FALSE)

cat("\nInterpretation guide:\n")
cat("  model_fewshot vs model_uncorrected -> did the k labels help the model?\n")
cat("  model_fewshot vs null_fewshot      -> FAIR test. Is the model worth\n")
cat("                                        anything once the null also gets\n")
cat("                                        the same k labels?\n")
cat("  offset_err                         -> how far off the k-sample offset\n")
cat("                                        estimate is (true offset SD = 4.83)\n")
