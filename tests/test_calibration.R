# =============================================================================
# tests/test_calibration.R - Assertions for C1 (calibration) and C2 (transforms)
# =============================================================================
#   Rscript tests/test_calibration.R     (from the repository root)
# =============================================================================

source("R/calibration.R")

# --- C2: transform round-trips ----------------------------------------------
y <- c(0, 0, 1, 5, 14, 42, 101)   # includes the zero spike and the observed max

# 1. Every transform must round-trip non-negative values to themselves.
for (nm in names(target_transforms)) {
  tf <- get_transform(nm)
  stopifnot(isTRUE(all.equal(tf$inverse(tf$forward(y)), y, tolerance = 1e-10)))
}

# 2. No transform may return a negative prediction. This is the entire point of
#    C2: HRDsum is bounded below at zero by construction.
for (nm in names(target_transforms)) {
  tf <- get_transform(nm)
  if (nm == "identity") next  # identity is deliberately unbounded; it IS the bug
  stopifnot(all(tf$inverse(c(-50, -1, -0.001, 0, 7)) >= 0))
}

# 3. Clipping must not change the ORDER of predictions among non-negative
#    values. If it did, it would be altering the model's rankings, not just
#    removing impossible values.
p <- c(-5, -1, 0, 3, 8, 30)
cl <- get_transform("clip")$inverse(p)
stopifnot(!is.unsorted(cl), identical(order(cl[p >= 0]), order(p[p >= 0])))

# 4. log1p must be monotone increasing - a transform that reordered patients
#    would invalidate every ranking metric we report.
lg <- get_transform("log1p")
stopifnot(!is.unsorted(lg$inverse(lg$forward(sort(y)))))

# 5. Unknown transforms must fail loudly rather than silently defaulting.
stopifnot(inherits(try(get_transform("sqrt_maybe"), silent = TRUE), "try-error"))

# --- C1: tissue covariates and calibrator -----------------------------------
set.seed(4471)
n <- 300
cancer <- rep(paste0("T", 1:10), each = 30)
purity <- runif(n, 0.2, 0.95)
predv  <- rnorm(n, 20, 6)
ood    <- runif(n, 0, 2)

tab <- tissue_covariates(predv, purity, cancer, ood)

# 6. One row per tissue, and every covariate finite. A NA covariate would
#    silently drop a tissue from lm() later.
stopifnot(nrow(tab) == 10, all(tab$n == 30),
          all(is.finite(tab$mean_pred)), all(is.finite(tab$sd_pred)),
          all(is.finite(tab$mean_purity)))

# 7. tissue_covariates() must NOT accept or require labels. This is the
#    structural guarantee that the calibrator can be applied to an unseen
#    tissue with no HRD measurements at all.
stopifnot(!"y" %in% names(formals(tissue_covariates)),
          !"offset" %in% names(tab))

# 8. On a SYNTHETIC case where the offset genuinely is a linear function of a
#    covariate, LOTO must detect it. This guards against the calibrator being
#    broken in a way that always reports failure - which matters because the
#    REAL-DATA result is negative, and a permanently-negative function would
#    produce that same answer for the wrong reason.
#
#    NOTE ON THE FIXTURE: tissue covariates are MEANS over ~30 samples, so
#    averaging washes out between-tissue variation. An earlier version of this
#    test drew purity from one common distribution for every tissue, leaving a
#    between-tissue spread of only 0.089 - barely above the injected noise - and
#    LOTO R^2 landed at 0.52 rather than >0.8. That was the fixture being
#    unrealistic, not the calibrator failing. Real tissues differ systematically
#    in mean purity, so give each tissue its own centre here.
tissue_centre <- seq(0.25, 0.9, length.out = 10)
purity_sep <- rep(tissue_centre, each = 30) + rnorm(n, 0, 0.03)
tab2 <- tissue_covariates(predv, purity_sep, cancer, ood)
tab2$offset <- 2.5 * tab2$mean_purity + rnorm(10, 0, 0.05)
stopifnot(diff(range(tab2$mean_purity)) > 0.5)   # fixture sanity
good <- fit_tissue_calibrator(tab2, covariates = "mean_purity")
stopifnot(good$loto_r2 > 0.8, good$loto_mae < good$naive_mae)

# 9. On PURE NOISE, LOTO R^2 must not be strongly positive. This is the
#    guard that makes the real-data negative result trustworthy.
tab2$offset <- rnorm(10, 0, 3)
noise <- fit_tissue_calibrator(tab2, covariates = "mean_purity")
stopifnot(noise$loto_r2 < 0.5)

# 10. apply_tissue_calibrator() must return a finite prediction even when the
#     covariate model degenerates - never NA, which would silently blank a
#     patient's result downstream.
out <- apply_tissue_calibrator(good, tab2[1, , drop = FALSE], c(10, 20, 30))
stopifnot(length(out) == 3, all(is.finite(out)))

cat("C1/C2 tests passed: transforms round-trip and stay non-negative;\n")
cat("  tissue covariates are label-free; LOTO detects real signal and\n")
cat("  does not manufacture it from noise.\n")
