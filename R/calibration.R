# =============================================================================
# R/calibration.R - Target transforms (C2) and cross-tissue calibration (C1).
# =============================================================================
#
# Addresses two defects measured in LOCO run 01 (docs/24):
#
#   C2  HRDsum is bounded at zero and 14.3% of labels are exactly 0, but the
#       elastic net is unbounded and emitted negative predictions for 3.2% of
#       samples (most negative: -11.53).
#
#   C1  The model mis-levels each tissue by 3.46 units on average (SD 4.83).
#       Removing that offset drops pooled MAE 9.049 -> 7.881, which is MORE
#       than the model's entire margin over the tissue-mean null (+0.085).
#
# These live outside R/model.R so that fit_en() keeps its current contract and
# the new behaviour can be unit tested and switched off independently.
#
# -----------------------------------------------------------------------------
# THE CENTRAL DIFFICULTY WITH C1 - READ BEFORE MODIFYING
# -----------------------------------------------------------------------------
# A per-tissue offset is normally estimated from labelled samples OF THAT
# TISSUE. In leave-one-cancer-out - and in the real pediatric application - the
# held-out tissue has NO LABELS BY DEFINITION. A lookup table keyed by cancer
# type is therefore useless exactly where it is needed, and worse, fitting one
# on the held-out tissue's own labels is leakage that invalidates the fold.
#
# The approach taken here: learn a REGRESSION from tissue-level COVARIATES to
# the tissue's offset, using the 29 training tissues as 29 observations. For an
# unseen tissue we compute the same covariates (all of which are available
# WITHOUT labels, from the beta matrix and metadata alone) and predict its
# offset. Whether this actually generalises is an empirical question, which is
# why fit_tissue_calibrator() reports leave-one-tissue-out CV error on the
# training tissues rather than in-sample fit.
#
# HONEST WARNING: with only 29 observations and a handful of covariates this is
# a small-n regression and is easy to overfit. The LOTO-CV figure is the number
# to trust. If LOTO R^2 is near zero, the correct conclusion is that the offset
# is NOT predictable from these covariates, and the fallback is to report
# within-tissue relative ranks instead of absolute values (docs/24 section 6).
# =============================================================================


# -----------------------------------------------------------------------------
# C2: target transforms
# -----------------------------------------------------------------------------
# Each transform is a matched (forward, inverse) pair plus a name. Keeping them
# paired in one object prevents the classic bug of transforming the target and
# forgetting to back-transform the prediction.
#
# "identity"  fit on the raw scale. The run 01 baseline.
# "clip"      fit on the raw scale, then clip predictions at 0. This is a
#             POST-HOC transform: it needs no refit, cannot change the ranking
#             of samples (Spearman with the raw prediction is exactly 1), and
#             only removes a priori impossible values.
# "log1p"     fit on log(1 + y), back-transform with expm1. This changes the
#             LOSS FUNCTION: absolute errors on high-HRD tumours are
#             down-weighted relative to low-HRD ones. That is a real scientific
#             trade-off, not a free improvement, because high-HRD tumours are
#             the clinically actionable ones. Evaluate before adopting.
#
# NOTE on the zero spike: log1p maps 0 -> 0, so the 14.3% of samples at exactly
# zero remain a point mass that a continuous regressor cannot reproduce. log1p
# addresses the RIGHT SKEW, not the zero inflation. Genuine zero inflation would
# need a hurdle/two-part model, which is out of scope here and recorded as
# future work in docs/24.
target_transforms <- list(
  identity = list(name = "identity",
                  forward = function(y) y,
                  inverse = function(p) p),
  clip     = list(name = "clip",
                  forward = function(y) y,
                  inverse = function(p) pmax(p, 0)),
  log1p    = list(name = "log1p",
                  forward = function(y) log1p(y),
                  # expm1 can still return a small negative if p < 0, so clip
                  # after back-transforming. Belt and braces.
                  inverse = function(p) pmax(expm1(p), 0))
)

get_transform <- function(name) {
  if (!name %in% names(target_transforms)) {
    stop(sprintf("Unknown transform '%s'. Available: %s",
                 name, paste(names(target_transforms), collapse = ", ")))
  }
  target_transforms[[name]]
}


# -----------------------------------------------------------------------------
# C1: tissue-level covariates
# -----------------------------------------------------------------------------
# Build one row per tissue describing that tissue WITHOUT using its HRD labels.
#
# Every covariate here must be computable for a tissue we have never seen and
# have no labels for. That constraint is what makes the calibrator transferable,
# and it is the reason this function deliberately does NOT accept `y`.
#
#   mean_pred     the model's own mean prediction for that tissue. Available at
#                 inference time because it needs no labels. This is the single
#                 most informative covariate: if the model systematically
#                 over-predicts a tissue, its mean prediction tends to sit high
#                 relative to the training distribution.
#   sd_pred       within-tissue spread of predictions. A tissue the model thinks
#                 is homogeneous behaves differently from a heterogeneous one.
#   mean_purity   tumour purity, from metadata, no labels needed.
#   sd_purity     purity heterogeneity.
#   mean_ood      mean out-of-distribution distance. Large values flag a tissue
#                 whose methylation sits far from the training manifold, which
#                 is precisely when we expect calibration to drift.
#   n             cohort size, as a precision weight.
tissue_covariates <- function(pred, purity, cancer, ood_score = NULL) {
  stopifnot(length(pred) == length(cancer))
  df <- data.frame(pred = pred, purity = purity, cancer = as.character(cancer),
                   ood = if (is.null(ood_score)) NA_real_ else ood_score,
                   stringsAsFactors = FALSE)
  sp <- split(df, df$cancer)
  out <- do.call(rbind, lapply(names(sp), function(k) {
    s <- sp[[k]]
    data.frame(
      cancer_type = k,
      n           = nrow(s),
      mean_pred   = mean(s$pred, na.rm = TRUE),
      sd_pred     = stats::sd(s$pred, na.rm = TRUE),
      mean_purity = mean(s$purity, na.rm = TRUE),
      sd_purity   = stats::sd(s$purity, na.rm = TRUE),
      mean_ood    = mean(s$ood, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  # A single-sample tissue has undefined SD. Substitute 0 rather than dropping
  # the tissue, so the caller keeps a complete table.
  out$sd_pred[!is.finite(out$sd_pred)] <- 0
  out$sd_purity[!is.finite(out$sd_purity)] <- 0
  rownames(out) <- NULL
  out
}


# -----------------------------------------------------------------------------
# fit_tissue_calibrator(): learn offset ~ tissue covariates across tissues.
# -----------------------------------------------------------------------------
# `tissue_tab` must contain one row per TRAINING tissue, with the columns from
# tissue_covariates() plus `offset` = mean(prediction - truth) for that tissue.
#
# The model is a plain linear regression. With ~29 observations anything more
# flexible would fit noise. Covariates are chosen by the caller; the default
# set is deliberately small for the same reason.
#
# Returns the fitted model plus a LEAVE-ONE-TISSUE-OUT cross-validated R^2 and
# MAE. The LOTO figures are the ones that matter: in-sample R^2 on 29 points
# with 4 predictors will look good even if the relationship is noise.
fit_tissue_calibrator <- function(tissue_tab,
                                  covariates = c("mean_pred","sd_pred","mean_purity","mean_ood")) {
  stopifnot("offset" %in% names(tissue_tab))
  # Drop covariates that are constant or all-NA - they cannot contribute and
  # would make lm() rank-deficient.
  usable <- covariates[vapply(covariates, function(cv) {
    v <- tissue_tab[[cv]]
    !is.null(v) && sum(is.finite(v)) > 2 && stats::sd(v[is.finite(v)]) > 0
  }, logical(1))]
  if (length(usable) == 0) stop("No usable tissue covariates")

  f <- stats::as.formula(paste("offset ~", paste(usable, collapse = " + ")))
  fit <- stats::lm(f, data = tissue_tab)

  # LEAVE-ONE-TISSUE-OUT CV. This mirrors the deployment situation exactly:
  # predict the offset of a tissue that contributed nothing to the calibrator.
  n <- nrow(tissue_tab)
  loto <- vapply(seq_len(n), function(i) {
    m <- try(stats::lm(f, data = tissue_tab[-i, , drop = FALSE]), silent = TRUE)
    if (inherits(m, "try-error")) return(NA_real_)
    as.numeric(stats::predict(m, newdata = tissue_tab[i, , drop = FALSE]))
  }, numeric(1))

  ok <- is.finite(loto) & is.finite(tissue_tab$offset)
  ss_res <- sum((tissue_tab$offset[ok] - loto[ok])^2)
  ss_tot <- sum((tissue_tab$offset[ok] - mean(tissue_tab$offset[ok]))^2)

  list(
    model      = fit,
    covariates = usable,
    formula    = f,
    # Baseline to beat: predicting every tissue's offset as the global mean
    # offset. If LOTO MAE is not below this, the covariates add nothing.
    loto_pred  = loto,
    loto_mae   = mean(abs(tissue_tab$offset[ok] - loto[ok])),
    naive_mae  = mean(abs(tissue_tab$offset[ok] - mean(tissue_tab$offset[ok]))),
    loto_r2    = if (ss_tot > 0) 1 - ss_res/ss_tot else NA_real_,
    n_tissues  = sum(ok)
  )
}


# apply_tissue_calibrator(): predict and subtract the offset for a new tissue.
# `new_tab` is a one-row tissue_covariates() table for the unseen tissue.
# Returns corrected predictions. Never call this with an offset derived from
# the held-out tissue's own labels - that is the leakage this design avoids.
apply_tissue_calibrator <- function(cal, new_tab, pred) {
  est <- as.numeric(stats::predict(cal$model, newdata = new_tab))
  if (!is.finite(est)) est <- 0  # fall back to no correction, never to NA
  pred - est
}
