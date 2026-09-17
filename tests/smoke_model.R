# =============================================================================
# tests/smoke_model.R - Fast correctness checks for the R modelling core.
# =============================================================================
#
# USAGE
#   Rscript tests/smoke_model.R          (from the repository root)
#
# This is a stopifnot() script, not a testthat suite, so it is NOT picked up by
# `python -m unittest discover -s tests`. Run it explicitly. It uses a tiny
# synthetic matrix and finishes in seconds, so there is no excuse for skipping
# it before a long training run.
#
# WHY EACH CHECK EXISTS - these are regression guards for specific, realistic
# failure modes, not filler assertions:
#
#   1. Fit-object immutability. apply_preprocess() must not mutate the transform
#      it is handed. If it ever did, applying a model to held-out data would
#      silently alter the model - the exact shape of a leakage bug.
#   2. Conformal degradation. With too few calibration points the interval must
#      come back Inf ("no usable interval") rather than a falsely narrow one.
#   3. Single-sample invariance. A prediction for one sample alone must equal
#      that sample's prediction within a batch. If this fails, some statistic is
#      being computed across the input batch at predict time - again, leakage.
#      This is the check that matters most for a single-sample clinical assay.
#   4. Column-order invariance. Feeding probes in reverse order must give an
#      identical answer, proving apply_preprocess() aligns features by NAME
#      rather than by position. A positional bug here would silently pair every
#      probe with the wrong coefficient.
#   5. Duplicate-patient rejection. Must be a hard error, since repeated
#      patients leak across fold boundaries.
#   6. Degenerate-metric handling. A constant predictor yields NA R2, not a
#      crash or a misleading number.
#
# NOT COVERED HERE: any of the leakage boundaries at the SCRIPT level (that
# locked_CNS never enters training, that the calibration set is disjoint from
# the training set), and nothing at all about the Python pipeline that builds
# beta.tsv. Passing this file does not mean the pipeline is correct.
# =============================================================================

source("R/model.R")

# --- Synthetic fixture -----------------------------------------------------
# 90 samples x 30 probes of uniform noise in [0,1], mimicking beta values.
# Fixed seed keeps the whole script deterministic.
set.seed(26);x<-matrix(runif(90*30),90,30,dimnames=list(paste0("s",1:90),paste0("cg",1:30)))
# y depends on only the first two probes plus noise, so a working elastic net
# has genuine signal to find while the other 28 probes act as distractors.
# 90 unique patients (no duplicates) across 3 cancer types, 30 samples each -
# enough for inner_folds() to take its primary cancer-grouped path.
y<-30*x[,1]-15*x[,2]+rnorm(90);patients<-paste0("p",1:90);cancer<-rep(c("A","B","C"),each=30)

# --- Check 1: the fit object is not mutated by being applied ---------------
# Learn the transform on the first 60 rows, snapshot it, apply it to rows 61-90,
# then confirm the snapshot still matches.
pp<-fit_preprocess(x[1:60,,drop=FALSE],10);before<-pp
z<-apply_preprocess(x[61:90,,drop=FALSE],pp)
stopifnot(identical(pp,before),identical(colnames(z),pp$features))

# --- Check 2: conformal interval degrades honestly -------------------------
# n=5 cannot support 95% coverage (ceiling((5+1)*0.95) = 6 > 5), so Inf.
# n=30 can (ceiling(31*0.95) = 30 <= 30), so a finite width.
stopifnot(is.infinite(conformal_q(1:5,1:5)),is.finite(conformal_q(1:30,1:30)))

# --- Fit a model for the remaining checks ----------------------------------
b<-fit_en(x,y,patients,cancer,max_features=20)

# --- Check 3: single-sample prediction == batch prediction ------------------
# Tolerance 1e-10 means bitwise-equivalent in practice. Any cross-sample
# statistic at predict time would break this immediately.
p1<-predict_en(b,x[1,,drop=FALSE]);pall<-predict_en(b,x)
stopifnot(isTRUE(all.equal(p1$predicted_reference_HRDsum,pall$predicted_reference_HRDsum[1],tolerance=1e-10)))

# --- Check 4: column order does not matter ---------------------------------
# Same sample, probes reversed. Must produce the identical prediction.
p2<-predict_en(b,x[1,rev(seq_len(ncol(x))),drop=FALSE])
stopifnot(isTRUE(all.equal(p1$predicted_reference_HRDsum,p2$predicted_reference_HRDsum,tolerance=1e-10)))

# --- Check 5: duplicate patients are rejected ------------------------------
# All 90 rows given the same patient ID must raise, via the anyDuplicated()
# guard in inner_folds().
try_duplicate<-try(fit_en(x,y,rep("same",90),cancer),silent=TRUE);stopifnot(inherits(try_duplicate,"try-error"))

# --- Check 6: R2 is NA when the reference values are constant ---------------
# sum((y - mean(y))^2) is 0, so R2 is undefined; metrics() must return NA
# rather than dividing by zero.
stopifnot(is.na(metrics(rep(1,5),1:5)$R2))

cat("R preprocessing/grouping/order/single-sample/conformal smoke tests passed\n")
