# =============================================================================
# scripts/train_baseline.R - Development training and leave-one-cancer-out
#                            evaluation for the methylation -> HRDsum model.
# =============================================================================
#
# USAGE
#   Rscript scripts/train_baseline.R <beta.tsv> <master_samples.tsv> <output_dir>
#
#   e.g. Rscript scripts/train_baseline.R \
#          data/processed/beta.tsv \
#          data/processed/master_samples.tsv \
#          results/baseline
#
#   Must be run from the repository root: source("R/model.R") below is a
#   relative path with no fallback resolution.
#
# WHAT THIS SCRIPT DOES, IN ORDER
#   1. Load and validate the beta matrix and the master sample table.
#   2. Refuse to run on an engineering/smoke fixture (provenance gate).
#   3. Hold out GBM and LGG entirely - the "locked CNS" partition.
#   4. Leave-one-cancer-out over the remaining development cancers, writing
#      per-fold predictions, metrics and fitted bundles.
#   5. Reserve a stratified 20% calibration set, fit a final model on the rest,
#      derive a split-conformal interval, and freeze that model to disk.
#
# WHAT IT DELIBERATELY DOES NOT DO
#   It never evaluates on the locked CNS partition. Architecture and
#   hyperparameter decisions must be finalised BEFORE anything touches the
#   locked set, otherwise the locked set stops being a held-out test and
#   becomes just another tuning signal. See the note at the bottom.
#
# ---------------------------------------------------------------------------
# OPERATIONAL WARNING - THIS SCRIPT WILL NOT RUN AS-IS AT FULL SCALE
# ---------------------------------------------------------------------------
# Measured against the real data/processed/beta.tsv (336,480 probes x 7,707
# samples):
#
#   MEMORY: line 6's fread() materialises ~21 GB as a data.frame, and the
#   transpose chain on the x <- t(as.matrix(...)) line creates three to four
#   live copies before modelling starts - 60 GB or more at peak. `beta` is
#   never rm()'d. Submit to a large-memory LSF queue; do not run on a login
#   node (typically ~3 GB free).
#
#   TIME: fit_preprocess() in R/model.R uses apply() over ~336k columns,
#   benchmarked at ~2.8 min per call. Nested LOCO invokes it roughly 900 times
#   (30 inner folds x 30 outer folds) -> about 42 hours. Substituting
#   matrixStats::colMedians/colVars brings that to ~17 hours.
#
# Recommended fixes before a real run: convert the matrix to a float HDF5 or
# bigstatsr FBM backing (halves it to ~10 GB), swap apply() for matrixStats,
# and add rm(beta); gc() after the transpose.
#
# ---------------------------------------------------------------------------
# AUDIT NOTES ON INTERPRETING THE OUTPUT (2026-09-16)
# ---------------------------------------------------------------------------
#   * The ONLY null comparator produced is null_MAE, the training-mean
#     prediction. There is no permuted-label control and no cancer-type-mean
#     (tissue identity) control. HRDsum varies strongly by tissue, so without a
#     tissue null you cannot tell whether the model learned HRD biology or
#     learned to recognise the tumour type. Add these before presenting results.
#   * null_MAE uses the TRAINING mean as its baseline, while the R2 reported by
#     metrics() uses the HELD-OUT CANCER's own mean. They are different
#     baselines; do not compare them side by side without saying so.
#   * `selected_features` counts features ENTERING glmnet, not features with
#     non-zero coefficients. It is not a sparsity measure.
#   * The LOCO metrics and the frozen model come from DIFFERENT fits - the
#     frozen model excludes the 20% calibration reservation. macro_metrics.txt
#     does not describe the artefact that ships in frozen_nonCNS.rds.
#   * purity and ploidy are present in the master table and never used here,
#     despite being plausible confounders (methylation predicts tumour purity;
#     purity correlates with scar scores).
# =============================================================================

# Run from repository root: Rscript scripts/train_baseline.R beta.tsv master.tsv results/run01
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("Usage: Rscript scripts/train_baseline.R beta.tsv master.tsv output_dir")
# All modelling functions live in R/model.R. Relative path => run from repo root.
source("R/model.R")
if (!requireNamespace("data.table",quietly=TRUE)) stop("Install dependencies with scripts/setup.R")

# --- Load the beta matrix --------------------------------------------------
# On-disk layout is probes in ROWS, samples in COLUMNS, first column = probe_id.
# check.names=FALSE preserves TCGA barcodes verbatim; R would otherwise mangle
# the hyphens in "TCGA-A2-A0SV-01A-..." into dots.
# MEMORY: this is the ~21 GB allocation described in the header.
beta <- data.table::fread(args[1],data.table=FALSE,check.names=FALSE)
if (anyDuplicated(beta[[1]])) stop("Duplicate probe IDs")

# Transpose to the samples-in-rows orientation that R/model.R expects, then
# label columns with probe IDs so apply_preprocess() can match features BY NAME.
# storage.mode ensures a numeric matrix even if fread typed a column as integer
# or character. NOTE: if any cell were non-numeric this coercion emits a warning
# ("NAs introduced by coercion") rather than an error, and those NAs would then
# flow onward looking like legitimate missing data.
x <- t(as.matrix(beta[,-1,drop=FALSE])); storage.mode(x)<-"double";colnames(x)<-beta[[1]]

# Beta values are methylation proportions and must lie in [0,1]. The is.finite()
# conjunct means NA is tolerated (handled later by imputation) but an out-of-range
# real number is fatal. Caveat: Inf passes this check, because is.finite(Inf) is
# FALSE - see the related note in R/model.R's fit_preprocess().
if(any(is.finite(x)&(x<0|x>1))) stop("Expected beta values in [0,1]")

# --- Load and validate the sample metadata ---------------------------------
meta <- read.delim(args[2],check.names=FALSE,stringsAsFactors=FALSE)
# One row per patient AND per specimen. Duplicates would leak a patient across
# the fold boundary during cross-validation.
if(anyDuplicated(meta$sample_id)||anyDuplicated(meta$patient_id)) stop("Duplicate patient/specimen")
# Every matrix column must have metadata. The reverse is allowed: extra metadata
# rows (samples not in this matrix) are simply dropped by the match() below.
if(any(!rownames(x)%in%meta$sample_id)) stop("Missing sample metadata")
# Align metadata row order to matrix row order. EVERY subsequent meta$... lookup
# depends on this alignment holding, so it must come before any modelling.
meta <- meta[match(rownames(x),meta$sample_id),]

# Target sanity: HRDsum is a count of genomic scars, so it cannot be negative.
if(any(!is.finite(meta$HRDsum))||any(meta$HRDsum<0)) stop("Invalid target")
# Label integrity: HRDsum must equal the sum of its three components. This is
# the single most valuable check in the file - it catches column mis-mapping
# between HRD_LOH / LST / TAI, which would otherwise be invisible and would
# silently train the model against a scrambled target.
if(!all(abs(meta$HRDsum-meta$HRD_LOH-meta$LST-meta$TAI)<1e-6)) stop("Label sum mismatch")

# WARNING - THIS GATE IS CURRENTLY TAUTOLOGICAL.
# scripts/build_master.py writes quality_annotation as a hardcoded constant
# string for every row, so this condition can never be false. It reads like a QC
# check but verifies nothing. Replace it with a real provenance assertion (for
# example, comparing a recorded source checksum) if you want an actual gate.
if (!all(meta$quality_annotation=="published_450K_no_exclusion")) stop("Unadjudicated QC")

# --- Provenance gate: refuse to train on a smoke fixture --------------------
# scripts/prepare_beta.py writes a sidecar <matrix_basename>.provenance.json
# recording whether the matrix was built with a real probe allowlist or with the
# tiny --smoke-probes engineering subset. Without this gate it would be
# alarmingly easy to present fixture output as a scientific result.
# Never silently run a smoke-probe fixture as scientific evidence.
prov <- sub("\\.[^.]+$",".provenance.json",args[1])
if(!file.exists(prov)) stop("Missing matrix provenance JSON")
if(!requireNamespace("jsonlite",quietly=TRUE)) stop("Install jsonlite")
pr <- jsonlite::fromJSON(prov)
if(isTRUE(pr$engineering_only)) stop("Engineering-only matrix. Build a real technical probe allowlist first.")
if(is.null(pr$probe_allowlist_sha256)) stop("Missing technical allowlist hash")

out <- args[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)

# --- Define the locked CNS partition ---------------------------------------
# GBM and LGG (642 samples) are withheld as the adult-CNS stand-in for the
# eventual pediatric high-grade glioma transfer question. They are never fitted
# and never scored in this script.
#
# NOTE: this re-derives the lock from a hardcoded cancer-type list rather than
# reading the `partition` column that build_master.py already wrote. The two
# agree today (verified: 642 locked_CNS rows, exactly GBM+LGG, zero anomalies),
# but they are two sources of truth that could drift. Reading partition and
# stopifnot()-ing agreement would be safer.
cns <- meta$cancer_type %in% c("GBM","LGG")

# Anti-fixture guards: refuse to make transfer claims from a toy cohort.
if(length(unique(meta$cancer_type[!cns]))<3) stop("Need >=3 development cancers; six-patient smoke subset is insufficient")
if(sum(!cns)<50) stop("Need >=50 development patients for this protocol; do not claim transfer from a tiny fixture")

# ===========================================================================
# PHASE 1: Leave-one-cancer-out (LOCO) development evaluation
# ===========================================================================
# Each iteration holds out one entire cancer type. This tests whether the model
# has learned something transferable about HRD biology, or has merely learned
# tissue-specific methylation signatures. Holding out random SAMPLES instead
# would let the model recognise the tissue and look far better than it is.
all_predictions <- list();all_metrics<-list()
for (type in sort(unique(meta$cancer_type[!cns]))) {
 # Both masks carry !cns, so locked CNS samples appear in neither train nor test.
 tr <- !cns & meta$cancer_type!=type;te <- !cns & meta$cancer_type==type

 # fit_en() runs its own inner CV for hyperparameters using only these training
 # rows - the held-out cancer is invisible to hyperparameter selection.
 b <- fit_en(x[tr,,drop=FALSE],meta$HRDsum[tr],meta$patient_id[tr],meta$cancer_type[tr])

 p <- predict_en(b,x[te,,drop=FALSE]);p$actual<-meta$HRDsum[te];p$cancer_type<-type;p$patient_id<-meta$patient_id[te]
 # The null comparator: predict this fold's training mean for every sample.
 # Stored per-row so downstream analysis can compute any null metric it wants.
 p$null_prediction<-b$training_mean;p$split<-"development_LOCO"
 all_predictions[[type]]<-p

 # Metrics are computed over ALL held-out samples, including those flagged
 # non-reportable. That is the non-inflating choice: reporting only "confident"
 # predictions would be selective reporting. abstention_rate is tracked
 # separately so the trade-off stays visible.
 mm<-metrics(p$actual,p$predicted_reference_HRDsum);mm$cancer_type<-type;mm$train_n<-sum(tr);mm$selected_features<-length(b$preprocess$features);mm$alpha<-b$alpha;mm$lambda<-b$lambda;mm$null_MAE<-mean(abs(p$actual-b$training_mean));mm$abstention_rate<-mean(!p$reportable)
 all_metrics[[type]]<-mm

 # Persist the per-fold bundle and inner-fold assignments so the grid search and
 # fold structure remain auditable after the run.
 saveRDS(b,file.path(out,paste0("loco_",type,".rds")))
 write.table(b$inner_folds,file.path(out,paste0("inner_folds_",type,".tsv")),sep="\t",row.names=FALSE,quote=FALSE)
}

write.table(do.call(rbind,all_predictions),file.path(out,"loco_predictions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
mt<-do.call(rbind,all_metrics);write.table(mt,file.path(out,"loco_metrics.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
# Macro average: every cancer counts equally regardless of sample count, so BRCA
# (n=743) cannot drown out CHOL (n=35). Propagates NA loudly if any fold failed.
writeLines(paste("macro_MAE",mean(mt$MAE)),file.path(out,"macro_metrics.txt"))

# ===========================================================================
# PHASE 2: Fit and freeze the final model, with a conformal interval
# ===========================================================================
# Deterministic calibration reservation within non-CNS cancers; no refit on calibration.
#
# Stratified 20% per cancer type, so the calibration set mirrors the development
# cohort's tissue composition rather than being dominated by the largest cancers.
# The seed makes the split reproducible.
#
# LATENT BUG: sample(ii, k) where ii has length 1 is interpreted by R as
# sample(1:ii, k) - it returns a random integer in 1..ii rather than ii itself.
# That stray index would point at an arbitrary row (possibly a locked_CNS one).
# It does not fire on the current cohort because the smallest development group
# is OV with n=10 (floor(10*0.2) = 2), but it will fire on any subset or
# re-partition that produces a singleton group. Fix with an explicit
# length(ii)==1L guard.
set.seed(260910);dev<-which(!cns);cal<-unlist(lapply(split(dev,meta$cancer_type[dev]),function(ii) sample(ii,max(1L,floor(length(ii)*0.2)))))
tr<-setdiff(dev,cal)

# Fit on development-minus-calibration only. The calibration samples must stay
# unseen or the conformal interval below loses its validity guarantee.
b<-fit_en(x[tr,,drop=FALSE],meta$HRDsum[tr],meta$patient_id[tr],meta$cancer_type[tr])

# Split-conformal: predict the untouched calibration set, then take the residual
# quantile as the interval half-width. No refit happens on `cal`.
cp<-predict_en(b,x[cal,,drop=FALSE]);b$interval_q<-conformal_q(meta$HRDsum[cal],cp$predicted_reference_HRDsum)

# Record exactly which patients were used for what, so any later question about
# contamination can be answered from the frozen object itself.
b$calibration_ids<-meta$patient_id[cal];b$training_ids<-meta$patient_id[tr];b$probe_provenance<-pr
b$interval_note<-"Split-conformal residual interval; no pediatric/domain-shift coverage guarantee"

# THE SHIPPING ARTEFACT. scripts/predict_frozen.R consumes this file.
saveRDS(b,file.path(out,"frozen_nonCNS.rds"))

# Export coefficients for interpretation. Rows with coefficient 0 were shrunk
# out by the elastic net; the non-zero count here is the real sparsity measure
# (unlike the `selected_features` column in loco_metrics.tsv).
cc<-as.matrix(coef(b$model,s=b$lambda))
write.table(data.frame(feature=rownames(cc),coefficient=as.numeric(cc[,1]),stringsAsFactors=FALSE),file.path(out,"frozen_coefficients.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Full role assignment for every sample: locked_CNS / calibration / training.
# This is the audit trail proving the locked set was never fitted.
write.table(data.frame(sample_id=meta$sample_id,patient_id=meta$patient_id,cancer_type=meta$cancer_type,role=ifelse(cns,"locked_CNS",ifelse(seq_len(nrow(meta))%in%cal,"calibration","training"))),file.path(out,"final_partitions.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

# Deliberately do not evaluate CNS here: architecture selection must precede locked testing.
# Scoring the locked set now would convert it from a held-out test into part of
# the development loop. Evaluate it exactly once, after the model is frozen, via
# scripts/predict_frozen.R.
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
cat("Development LOCO complete; frozen_nonCNS.rds saved. Review development and lock before external prediction.\n")
