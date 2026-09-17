#!/usr/bin/env Rscript
# =============================================================================
# scripts/benchmark_preprocess.R
# =============================================================================
#
# PURPOSE
#   Measure how long fit_preprocess() actually takes on the REAL beta matrix at
#   its full width, so that the projected wall time for a full nested
#   leave-one-cancer-out (LOCO) run rests on a measurement rather than an
#   extrapolation from small test matrices.
#
#   Earlier estimates were extrapolated from 20k- and 60k-probe benchmarks run
#   on a login node. Extrapolation assumes cost grows linearly with the number
#   of probes. That assumption is reasonable for compiled column reductions but
#   can break down once a matrix no longer fits comfortably in cache or RAM,
#   where memory bandwidth rather than arithmetic becomes the limit. This script
#   exists to check whether it holds.
#
# USAGE
#   Rscript scripts/benchmark_preprocess.R <beta.tsv> <master_samples.tsv> <out_dir>
#
#   Intended to be run through bsub (see scripts/lsf_benchmark_preprocess.bsub)
#   because loading the matrix needs far more memory than a login node provides.
#
# WHAT IT DOES NOT DO
#   No model is fitted and no result is produced. This is a timing harness only.
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop("Usage: benchmark_preprocess.R <beta.tsv> <master_samples.tsv> <out_dir>")

beta_path <- args[1]; meta_path <- args[2]; out_dir <- args[3]
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Report the machine we landed on. Timings are only interpretable alongside the
# hardware that produced them, and an LSF array can place jobs on different
# node types within the same queue.
cat("host          :", Sys.info()[["nodename"]], "\n")
cat("R version     :", as.character(getRversion()), "\n")
cat("matrixStats   :", if (requireNamespace("matrixStats", quietly=TRUE))
      as.character(packageVersion("matrixStats")) else "NOT INSTALLED", "\n")
cat("started       :", format(Sys.time()), "\n\n")

source("R/model.R")

# --- Load the matrix --------------------------------------------------------
# On disk the layout is probes in ROWS and samples in COLUMNS. R/model.R expects
# the opposite (samples in rows), so we transpose after reading. Both the read
# and the transpose are timed separately, because the transpose alone briefly
# doubles peak memory and is a real part of the job's cost.
cat("reading matrix (this is the ~21 GB allocation) ...\n")
t_read <- system.time({
  beta <- data.table::fread(beta_path, data.table=FALSE, check.names=FALSE)
})[["elapsed"]]
cat(sprintf("  fread: %.1f s  ->  %d probes x %d samples\n",
            t_read, nrow(beta), ncol(beta)-1L))

cat("transposing to samples-in-rows ...\n")
t_tr <- system.time({
  x <- t(as.matrix(beta[,-1,drop=FALSE]))
  storage.mode(x) <- "double"
  colnames(x) <- beta[[1]]
})[["elapsed"]]
# Free the on-disk-orientation copy immediately; keeping both costs ~21 GB extra.
rm(beta); invisible(gc(verbose=FALSE))
cat(sprintf("  transpose: %.1f s  ->  %d samples x %d probes\n", t_tr, nrow(x), ncol(x)))
cat(sprintf("  matrix in memory: %.1f GB\n\n", as.numeric(object.size(x))/1024^3))

# --- Reproduce a realistic training-set size --------------------------------
# A LOCO fold trains on every development sample except one cancer type, and an
# inner fold trains on a little less again. Timing the full cohort would
# overstate the per-call cost, so we subset to a representative training size.
meta <- read.delim(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
meta <- meta[match(rownames(x), meta$sample_id),]
is_cns <- meta$cancer_type %in% c("GBM","LGG")

# Representative inner-fold training set: development samples minus two cancer
# types (one held out by the outer loop, one by the inner loop).
set.seed(260910)
types <- sort(unique(meta$cancer_type[!is_cns]))
drop2 <- types[1:2]
train_rows <- which(!is_cns & !meta$cancer_type %in% drop2)
cat(sprintf("representative inner-fold training set: %d samples (%d dev types minus 2)\n\n",
            length(train_rows), length(types)))

# --- Time fit_preprocess ----------------------------------------------------
# Three repeats: the first can be inflated by page faults as pages are touched
# for the first time, so we report the median as the number to plan with.
reps <- 3L
times <- numeric(reps)
for (i in seq_len(reps)) {
  cat(sprintf("fit_preprocess replicate %d/%d ...\n", i, reps))
  times[i] <- system.time({
    pp <- fit_preprocess(x[train_rows,,drop=FALSE], max_features=5000L)
  })[["elapsed"]]
  cat(sprintf("  %.1f s  (kept %d features)\n", times[i], length(pp$features)))
}

# --- Time apply_preprocess --------------------------------------------------
# apply_preprocess() runs twice per inner fold (training and validation), so it
# is a meaningful share of the total and must be measured, not assumed.
cat("\napply_preprocess (train split) ...\n")
t_apply <- system.time(invisible(apply_preprocess(x[train_rows,,drop=FALSE], pp, FALSE)))[["elapsed"]]
cat(sprintf("  %.1f s\n", t_apply))

# --- Project the full Phase 1 LOCO wall time --------------------------------
# Call counts follow the structure of fit_en() in R/model.R:
#   - one fit_preprocess per inner fold, plus one final refit per outer fold
#   - apply_preprocess twice per inner fold
#   - glmnet once per (inner fold x alpha), plus one final refit
n_types    <- length(types)
n_inner    <- n_types - 1L
pp_per_fold  <- n_inner + 1L
app_per_fold <- 2L * n_inner
t_pp  <- median(times)
sec_per_fold <- pp_per_fold*t_pp + app_per_fold*t_apply
total_sec    <- n_types * sec_per_fold

cat("\n================ PROJECTION ================\n")
cat(sprintf("development cancer types        : %d\n", n_types))
cat(sprintf("fit_preprocess median           : %.1f s\n", t_pp))
cat(sprintf("apply_preprocess                : %.1f s\n", t_apply))
cat(sprintf("calls per outer fold            : %d fit_preprocess, %d apply_preprocess\n",
            pp_per_fold, app_per_fold))
cat(sprintf("per outer fold (excl. glmnet)   : %.1f min\n", sec_per_fold/60))
cat(sprintf("SERIAL total, %d folds           : %.1f h\n", n_types, total_sec/3600))
cat(sprintf("PARALLEL (30-task array)        : %.1f min + glmnet + I/O\n", sec_per_fold/60))
cat("============================================\n")

# Persist the measurements so the projection can be audited later rather than
# recovered from a log file that may be rotated away.
write.table(
  data.frame(host=Sys.info()[["nodename"]],
             n_samples_total=nrow(x), n_probes=ncol(x),
             n_train_rows=length(train_rows),
             fread_sec=t_read, transpose_sec=t_tr,
             fit_preprocess_rep1=times[1], fit_preprocess_rep2=times[2],
             fit_preprocess_rep3=times[3], fit_preprocess_median=t_pp,
             apply_preprocess_sec=t_apply,
             projected_per_fold_min=sec_per_fold/60,
             projected_serial_hours=total_sec/3600),
  file.path(out_dir,"preprocess_benchmark.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat("\nfinished      :", format(Sys.time()), "\n")
cat("peak RSS is reported by LSF in the job summary below.\n")
