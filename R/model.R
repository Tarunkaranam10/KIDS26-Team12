# =============================================================================
# R/model.R - Fold-safe elastic-net helpers for predicting reference HRDsum
#             from DNA methylation beta values.
# =============================================================================
#
# WHAT THIS FILE IS
# -----------------
# The modelling core of the project. It is deliberately a library of plain
# functions with no side effects: nothing here reads a file, writes a file, or
# looks at the command line. Callers (scripts/train_baseline.R and
# scripts/predict_frozen.R) own all I/O. That separation is what makes the
# leakage discipline below auditable.
#
# THE CENTRAL INVARIANT - READ THIS BEFORE EDITING ANYTHING
# ---------------------------------------------------------
# Every quantity that is *learned from data* must be estimated from training
# rows only, stored inside a fit object, and then *applied* unchanged to
# held-out rows. In this file the learned quantities are:
#
#   1. which probes survive the missingness filter   (fit_preprocess: keep)
#   2. the per-probe median used for imputation      (fit_preprocess: med)
#   3. which probes survive variance ranking         (fit_preprocess: ii)
#   4. the per-probe centre and scale                (fit_preprocess: mu, s)
#   5. the elastic-net coefficients                  (fit_en: mod)
#   6. the out-of-distribution distance cutoff       (fit_en: ood_cut)
#
# All six live in the list returned by fit_preprocess()/fit_en(). apply_preprocess()
# and predict_en() only ever *consume* them. If you find yourself computing a
# mean, median, variance, or quantile inside apply_preprocess() or predict_en(),
# you have introduced test-set leakage and the reported performance becomes
# meaningless. Do not do it.
#
# A subtle corollary: glmnet is always called with standardize=FALSE. glmnet's
# default is standardize=TRUE, which would re-standardise using whatever matrix
# it was handed - including validation rows - silently undoing invariant (4).
# The FALSE is load-bearing, not stylistic.
#
# WHAT IS BEING PREDICTED, AND WHAT THAT DOES NOT MEAN
# -----------------------------------------------------
# The target is HRDsum = HRD_LOH + LST + TAI, a genomic-scar burden score
# derived from allele-aware data (SNP6 arrays). Methylation arrays cannot
# observe loss of heterozygosity or telomeric allelic imbalance directly, so
# this model is a *predictor of an independently measured reference score*, not
# a re-derivation of it. Nothing produced here is a clinical HRD assay.
#
# KNOWN LIMITATIONS OF THIS IMPLEMENTATION (audited 2026-09-16)
# --------------------------------------------------------------
#   * fit_preprocess() uses apply(z, 2, ...) over ~336k columns. Benchmarked at
#     roughly 2.8 minutes per call at full scale, and train_baseline.R calls it
#     ~900 times across nested LOCO -> about 42 hours. Swapping in
#     matrixStats::colMedians/colVars cuts this to ~17 hours. This is the single
#     biggest performance blocker in the project.
#   * The lambda grid is five hardcoded values spanning four decades and is not
#     anchored to glmnet's data-derived lambda.max. There is no check that the
#     selected lambda is interior to the grid, so a boundary optimum would be
#     accepted silently.
#   * The only null comparator produced downstream is the training mean. There
#     is no permuted-label control and no tissue-identity (cancer-type-mean)
#     control anywhere in the repository. Until those exist, a good-looking
#     result cannot be distinguished from the model having learned tissue
#     lineage rather than HRD biology.
# =============================================================================


# -----------------------------------------------------------------------------
# fit_preprocess(): LEARN the preprocessing transform from training rows only.
# -----------------------------------------------------------------------------
# Input : x            numeric matrix, samples in ROWS, probes in COLUMNS,
#                      with colnames set to probe IDs (e.g. "cg00000029").
#         max_features how many probes to keep after variance ranking.
#         max_missing  a probe is dropped if more than this fraction of training
#                      samples have a non-finite value for it.
#
# Output: a list describing the transform. It contains NO sample data - only
#         per-probe summary statistics - so it is safe to save and ship.
#
# This function must only ever be handed TRAINING rows. Every caller below
# subsets x before calling it.
fit_preprocess <- function(x, max_features=5000L, max_missing=0.05) {
  # Guard the shape contract up front. Duplicate probe names would make the
  # name-based column reconstruction in apply_preprocess() ambiguous, so reject
  # them here rather than producing a silently mis-aligned matrix later.
  stopifnot(is.matrix(x), !is.null(colnames(x)), !anyDuplicated(colnames(x)))

  # LEARNED QUANTITY 1: the missingness filter.
  # !is.finite() catches NA, NaN, Inf and -Inf together. colMeans() of that
  # logical matrix gives the per-probe missing fraction across TRAINING samples.
  keep <- colMeans(!is.finite(x)) <= max_missing
  if (!any(keep)) stop("No features pass train-only missingness")
  z <- x[,keep,drop=FALSE]

  # LEARNED QUANTITY 2: per-probe imputation medians, from training rows only.
  # Median rather than mean because beta values are bounded in [0,1] and often
  # strongly bimodal (methylated vs unmethylated), where a mean lands in a
  # trough that no real sample occupies.
  #
  # CAVEAT: na.rm=TRUE strips NA and NaN but NOT Inf. A column containing Inf
  # is therefore treated as missing by the mask above while still contributing
  # Inf to its own median. In practice such a column is removed by the v>0 &
  # is.finite(v) variance filter below, but the inconsistency is real.
  med <- apply(z,2,median,na.rm=TRUE)
  # Impute in place, column by column, using the medians just learned.
  for (j in seq_len(ncol(z))) z[!is.finite(z[,j]),j] <- med[j]

  # LEARNED QUANTITY 3: unsupervised feature selection by variance.
  # This ranks probes by how much they vary across TRAINING samples. It never
  # looks at y, so it is not supervised selection and does not leak the label.
  # Constant probes (v == 0) are dropped because they carry no information and
  # would produce a divide-by-zero at the scaling step below.
  v <- apply(z,2,var); ii <- which(is.finite(v) & v>0)
  # Sort by decreasing variance, breaking ties on probe NAME. The name tie-break
  # makes selection fully deterministic: two probes with identical variance
  # always resolve the same way regardless of column order in the input file.
  ii <- ii[order(-v[ii],names(v)[ii])]
  ii <- head(ii,max_features)
  if (length(ii)<2) stop("Fewer than two nonconstant features")

  # LEARNED QUANTITIES 4: centre and scale, again from training rows only.
  # Because the v>0 filter already ran, every retained column has sd > 0 and the
  # division in apply_preprocess() is safe.
  z <- z[,ii,drop=FALSE]; mu <- colMeans(z); s <- apply(z,2,sd)

  # Return the transform. med is subset to the finally-selected features so the
  # stored vectors are all the same length and in the same order as $features.
  list(features=colnames(z), median=med[colnames(z)],center=mu,scale=s,
       max_missing=max_missing)
}


# -----------------------------------------------------------------------------
# apply_preprocess(): APPLY a previously learned transform. Learns nothing.
# -----------------------------------------------------------------------------
# This is the held-out side of the leakage boundary. It reads pp$median,
# pp$center and pp$scale and never recomputes them.
#
# reject_missing controls whether excessive missingness is a hard error (used by
# the smoke test) or a soft flag (used everywhere in production, where
# predict_en() reports qc_fail instead of aborting the whole run for one bad
# sample).
apply_preprocess <- function(x,pp, reject_missing=TRUE) {
  if (anyDuplicated(colnames(x))) stop("Duplicate input features")

  # Rebuild the feature matrix BY NAME, not by position. This is what makes the
  # model invariant to column order in the incoming file, and it is explicitly
  # tested in tests/smoke_model.R by feeding in a reversed matrix.
  # Probes the model expects but the new data lacks stay NA here and are filled
  # with the stored training median below - which is the correct behaviour for a
  # single-sample assay that may not cover every probe.
  z <- matrix(NA_real_,nrow(x),length(pp$features),dimnames=list(rownames(x),pp$features))
  common <- intersect(colnames(x),pp$features);z[,common] <- x[,common,drop=FALSE]

  # Per-SAMPLE missingness (rowMeans, not colMeans): "how much of what this
  # model needs is absent for this particular sample?" This drives the
  # per-sample QC flag rather than any feature-level decision.
  miss <- rowMeans(!is.finite(z))
  if (reject_missing && any(miss>pp$max_missing)) stop("Too many missing model features")

  # Impute with the TRAINING medians. pp$median[j] is positionally aligned with
  # pp$features[j] because fit_preprocess() subset it by colnames(z).
  for (j in seq_len(ncol(z))) z[!is.finite(z[,j]),j] <- pp$median[j]

  # Centre and scale with the TRAINING statistics. Note this is the only
  # standardisation that ever happens - glmnet is called with standardize=FALSE.
  z <- sweep(sweep(z,2,pp$center,"-"),2,pp$scale,"/")

  # Carry the missing fraction along as an attribute so predict_en() can flag
  # samples without needing to recompute it.
  attr(z,"missing_fraction") <- miss;z
}


# -----------------------------------------------------------------------------
# inner_folds(): assign INNER cross-validation folds for hyperparameter tuning.
# -----------------------------------------------------------------------------
# Called from inside fit_en(), which is itself called once per OUTER
# leave-one-cancer-out fold. So this produces the inner layer of a properly
# nested CV: hyperparameters are chosen without ever touching the outer held-out
# cancer.
inner_folds <- function(patient,cancer,seed=260910L) {
  # One specimen per patient is a hard requirement. Two samples from the same
  # patient in different folds would leak that patient's methylation profile
  # across the fold boundary.
  if (anyDuplicated(patient)) stop("Use one specimen per patient for this MVP")

  # PRIMARY PATH: group by cancer type, so each remaining cancer becomes its own
  # inner fold. Combined with the outer loop this is leave-one-cancer-out nested
  # inside leave-one-cancer-out. Hyperparameters are therefore selected for their
  # ability to generalise ACROSS TISSUES, which is the actual scientific question,
  # rather than for within-tissue fit.
  if (length(unique(cancer))>=3) return(match(cancer,sort(unique(cancer))))

  # FALLBACK PATH: random 3-fold over patients. Unreachable from
  # train_baseline.R, which already refuses to run with fewer than 3 development
  # cancers; it exists for the smoke test and for future small-cohort use.
  if (length(patient)<9) stop("Need >=9 training patients for patient-fold fallback")
  set.seed(seed); as.integer(sample(rep(1:3,length.out=length(patient))))
}


# -----------------------------------------------------------------------------
# fit_en(): tune and fit the elastic net. The heart of the modelling code.
# -----------------------------------------------------------------------------
# Two phases:
#   Phase 1 - inner CV over an (alpha, lambda) grid to pick hyperparameters.
#   Phase 2 - refit preprocessing and model on ALL supplied rows using the
#             winning hyperparameters.
#
# Everything passed in here is training data for the current outer fold. The
# outer held-out cancer is never visible to this function.
fit_en <- function(x,y,patient,cancer,max_features=5000L,seed=260910L) {
  if (!requireNamespace("glmnet",quietly=TRUE)) stop("Install glmnet via scripts/setup.R")
  stopifnot(length(y)==nrow(x),all(is.finite(y)),length(patient)==length(y))
  # A constant target makes R2 undefined and the fit meaningless.
  if (sd(y)==0) stop("Constant training target")

  folds <- inner_folds(patient,cancer,seed)

  # The hyperparameter grid.
  #   alpha  = elastic-net mixing. 1.0 is pure lasso (sparse), 0.1 is nearly
  #            ridge (dense, handles correlated probes better). Methylation
  #            probes are heavily correlated in blocks, so the low-alpha end
  #            matters here.
  #   lambda = regularisation strength, four decades wide.
  #
  # LIMITATION: these are fixed absolute values, not anchored to the data-derived
  # lambda.max that glmnet would compute. On standardised predictors with y in
  # roughly 0-75, lambda=100 is almost certainly an intercept-only model and
  # lambda=0.01 is near-OLS. Nothing below verifies that the winner is INTERIOR
  # to the grid, so a boundary optimum is accepted without complaint.
  grid <- expand.grid(alpha=c(0.1,0.5,1),lambda=c(0.01,0.1,1,10,100))
  losses <- matrix(NA_real_,nrow(grid),length(unique(folds)))

  # ---- Phase 1: inner cross-validation -------------------------------------
  for (f in sort(unique(folds))) {
    tr <- folds!=f;va <- !tr

    # CRITICAL: preprocessing is re-learned from scratch on THIS inner fold's
    # training rows. It is not hoisted out of the loop. Hoisting it would let
    # every inner validation fold see statistics computed from itself.
    pp <- fit_preprocess(x[tr,,drop=FALSE],max_features)
    ztr <- apply_preprocess(x[tr,,drop=FALSE],pp,FALSE)
    zva <- apply_preprocess(x[va,,drop=FALSE],pp,FALSE)

    for (a in unique(grid$alpha)) {
      # Fit the whole lambda path in one call - glmnet is far more efficient
      # doing this than being called once per lambda. Decreasing order is
      # glmnet's expected convention for warm starts.
      # standardize=FALSE: see the header note. ztr is already standardised
      # using training-only statistics and glmnet must not redo it.
      mod <- glmnet::glmnet(ztr,y[tr],alpha=a,lambda=sort(unique(grid$lambda),decreasing=TRUE),standardize=FALSE)
      ids <- which(grid$alpha==a)
      for (g in ids) {
        # s= selects one lambda from the fitted path.
        pred <- as.numeric(predict(mod,zva,s=grid$lambda[g]))
        # Equal weight to each inner-validation cancer, independent of size.
        # NOTE: on the primary path each inner fold IS a single cancer, so this
        # tapply collapses to a plain mean within the fold; the macro-averaging
        # actually comes from rowMeans(losses) below. The tapply only does real
        # work on the patient-fold fallback path.
        losses[g,match(f,sort(unique(folds)))] <- mean(tapply(abs(pred-y[va]),cancer[va],mean))
      }
    }
  }

  # Average each grid point's loss across inner folds and take the winner.
  # Because each fold is a cancer, this is a macro average over tissues: a
  # hyperparameter setting that is excellent on BRCA and terrible everywhere
  # else will lose to one that is uniformly decent.
  grid$inner_macro_mae <- rowMeans(losses)
  best <- which.min(grid$inner_macro_mae); pp <- fit_preprocess(x,max_features)

  # ---- Phase 2: refit on all supplied training rows -------------------------
  # Note that fit_preprocess() is called again on the FULL x above - the final
  # model gets preprocessing learned from all of its own training data, which is
  # correct and still excludes the outer held-out cancer.
  z <- apply_preprocess(x,pp,FALSE)
  mod <- glmnet::glmnet(z,y,alpha=grid$alpha[best],lambda=sort(unique(grid$lambda),decreasing=TRUE),standardize=FALSE)

  # LEARNED QUANTITY 6: a crude out-of-distribution cutoff.
  # dist is the root-mean-square standardised distance of each training sample
  # from the training centroid. Samples far outside the training cloud get
  # flagged at prediction time.
  # Heuristic score only: not a calibrated domain-shift test.
  dist <- sqrt(rowMeans(z^2));ood_cut <- as.numeric(quantile(dist,0.99,names=FALSE))

  # The returned bundle is self-contained: it carries the transform, the model,
  # the winning hyperparameters, the full tuning table (so the grid search is
  # auditable after the fact), the fold assignments, the OOD cutoff, the seed,
  # and the training mean that serves as the null comparator downstream.
  list(preprocess=pp,model=mod,alpha=grid$alpha[best],lambda=grid$lambda[best],
       tuning=grid,inner_folds=data.frame(patient_id=patient,cancer_type=cancer,fold=folds),
       ood_cut=ood_cut,seed=seed,training_n=length(y),training_mean=mean(y),target="reference_HRDsum")
}


# -----------------------------------------------------------------------------
# predict_en(): apply a fitted bundle to new samples. Learns nothing.
# -----------------------------------------------------------------------------
# Returns one row per input sample with the prediction plus QC flags. Note it
# returns the raw prediction even when QC fails; it is the CALLER's job to
# decide whether to display it (predict_frozen.R blanks the display column but
# deliberately retains the raw value for analysis).
predict_en <- function(bundle,x) {
  # FALSE here means "flag, do not abort" - one bad sample in a batch should not
  # kill predictions for the rest.
  z <- apply_preprocess(x,bundle$preprocess,FALSE)
  raw <- as.numeric(predict(bundle$model,z,s=bundle$lambda))

  # Same distance metric as used to derive ood_cut at fit time, applied with the
  # stored cutoff rather than a freshly computed one.
  score <- sqrt(rowMeans(z^2));missing <- attr(z,"missing_fraction")
  fail <- missing>bundle$preprocess$max_missing
  ood <- score>bundle$ood_cut

  # reportable is the conjunction: a sample must both pass QC and sit inside the
  # training distribution before its prediction should be shown to anyone.
  data.frame(sample_id=rownames(x),predicted_reference_HRDsum=raw,
             missing_fraction=missing,ood_score=score,ood=ood,qc_fail=fail,
             reportable=!(fail|ood),stringsAsFactors=FALSE)
}


# -----------------------------------------------------------------------------
# conformal_q(): split-conformal prediction interval half-width.
# -----------------------------------------------------------------------------
# Given labels and predictions on a CALIBRATION set that the model was NOT fit
# on, return the residual quantile q such that [pred - q, pred + q] has the
# requested marginal coverage.
#
# k = ceiling((n+1) * coverage) is the textbook finite-sample-valid form (the
# +1 accounts for the new test point). If k > n there are too few calibration
# points to certify the requested coverage, and returning Inf is the honest
# answer: "no usable interval" rather than a falsely narrow one.
#
# IMPORTANT: this guarantee is marginal and assumes the test point is
# exchangeable with the calibration set. Applying a TCGA-calibrated interval to
# pediatric samples breaks exchangeability, so coverage there is not guaranteed.
conformal_q <- function(y,pred,coverage=0.95) {
  good <- is.finite(y)&is.finite(pred);res <- sort(abs(y[good]-pred[good]));n <- length(res)
  k <- ceiling((n+1)*coverage)
  if (!n || k>n) return(Inf)
  res[k]
}


# -----------------------------------------------------------------------------
# metrics(): evaluation panel for one set of predictions.
# -----------------------------------------------------------------------------
# Non-finite pairs are dropped, but n is reported so any such drop is visible.
metrics <- function(y,p) {
  ok <- is.finite(y)&is.finite(p);y<-y[ok];p<-p[ok];n<-length(y)
  if (!n) return(data.frame(n=0,MAE=NA,RMSE=NA,R2=NA,Pearson=NA,Spearman=NA,calibration_intercept=NA,calibration_slope=NA,bias=NA))

  # Correlations and the calibration regression need genuine variation in both
  # vectors; guard so a constant predictor yields NA rather than an error.
  variable <- n>2 && sd(y)>0 && sd(p)>0

  # Calibration: regress OBSERVED on PREDICTED. A perfectly calibrated model has
  # intercept 0 and slope 1. Slope < 1 means predictions are over-dispersed
  # (too extreme); slope > 1 means they are shrunk toward the mean.
  cal <- if (variable) coef(lm(y~p)) else c(NA,NA)

  data.frame(n=n,MAE=mean(abs(y-p)),RMSE=sqrt(mean((y-p)^2)),
             # R2 IS COMPUTED AGAINST THE MEAN OF y AS PASSED IN. Callers pass
             # one held-out cancer at a time, so mean(y) is that cancer's OWN
             # mean. This is the strict, conservative baseline: scoring above
             # zero requires beating an oracle that already knows the held-out
             # cancer's average HRDsum. It is NOT the training mean, and it is
             # therefore a harder bar than the null_MAE reported by
             # train_baseline.R - be careful not to compare the two directly.
             R2=if(n>1 && sum((y-mean(y))^2)>0) 1-sum((y-p)^2)/sum((y-mean(y))^2) else NA,
             # Pearson measures linear agreement; Spearman measures rank
             # agreement and is the more robust of the two if the relationship
             # turns out to be monotone but not linear.
             Pearson=if(variable) cor(y,p) else NA,Spearman=if(variable) cor(y,p,method="spearman") else NA,
             # bias is signed mean error: positive means systematic over-prediction.
             calibration_intercept=unname(cal[1]),calibration_slope=unname(cal[2]),bias=mean(p-y))
}
