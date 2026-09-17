#!/usr/bin/env Rscript
# =============================================================================
# scripts/loco_merge.R - combine per-fold outputs into cohort-level results.
# =============================================================================
#
# PURPOSE
#   scripts/loco_one_fold.R writes one set of files per outer fold. This script
#   gathers them and produces exactly the files scripts/train_baseline.R would
#   have written from its serial loop, plus the pooled confounding controls.
#
#   Splitting merge from fit matters because the pooled analyses are the ones
#   that actually detect tissue and purity confounding, and they cannot run
#   until every fold is finished.
#
# USAGE
#   Rscript scripts/loco_merge.R <out_dir> <master.tsv>
#
# IT WILL REFUSE TO RUN if any fold is missing. A merged result computed from 27
# of 30 folds would silently be a different (and better-looking) analysis than
# the one specified, because failed folds are rarely failing at random.
# =============================================================================

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Usage: Rscript scripts/loco_merge.R <out_dir> <master.tsv>")
out_dir <- args[1]; meta_path <- args[2]
fold_dir <- file.path(out_dir, "folds")
if (!dir.exists(fold_dir)) stop("No folds directory at ", fold_dir)

source("R/model.R")

meta <- read.delim(meta_path, check.names=FALSE, stringsAsFactors=FALSE)
cns <- meta$cancer_type %in% c("GBM","LGG")
types <- sort(unique(meta$cancer_type[!cns]))

# --- Completeness check -----------------------------------------------------
present <- sub("^predictions_","",sub("\\.tsv$","",basename(Sys.glob(file.path(fold_dir,"predictions_*.tsv")))))
missing <- setdiff(types, present)
if (length(missing)) {
  stop("Missing folds: ", paste(missing, collapse=", "),
       "\nRerun those array indices before merging. Refusing to merge a partial cohort.")
}
cat("all", length(types), "folds present\n\n")

read_all <- function(pattern) {
  do.call(rbind, lapply(Sys.glob(file.path(fold_dir, pattern)),
                        function(f) read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)))
}
pooled    <- read_all("predictions_*.tsv")
mt        <- read_all("metrics_*.tsv")
all_nulls <- read_all("nullpanel_*.tsv")

write.table(pooled, file.path(out_dir,"loco_predictions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(mt, file.path(out_dir,"loco_metrics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(all_nulls, file.path(out_dir,"loco_null_panel_by_cancer.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
writeLines(paste("macro_MAE", mean(mt$MAE)), file.path(out_dir,"macro_metrics.txt"))

# Surface any fold whose tuning hit a path endpoint. One or two is unremarkable;
# a majority means the lambda path needs widening before the result is trusted.
if ("lambda_at_boundary" %in% names(mt)) {
  nb <- sum(as.logical(mt$lambda_at_boundary), na.rm=TRUE)
  cat(sprintf("lambda at path boundary in %d of %d folds\n", nb, nrow(mt)))
  if (nb > nrow(mt)/2) cat("  WARNING: most folds hit a boundary; widen the lambda path.\n")
}

# ===========================================================================
# POOLED CONFOUNDING CONTROLS
# ===========================================================================
# Per-fold, every sample shares one cancer type, so the tissue-mean null is just
# that fold's mean. Pooling restores between-tissue variation, which is the only
# setting where the tissue confound is visible. See docs/22.
pooled_null <- null_panel(pooled$actual, pooled$predicted_reference_HRDsum,
                          pooled$cancer_type, mean(meta$HRDsum[!cns]))
write.table(pooled_null, file.path(out_dir,"pooled_null_panel.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

pooled_within <- within_tissue_metrics(pooled$actual, pooled$predicted_reference_HRDsum, pooled$cancer_type)
write.table(pooled_within, file.path(out_dir,"pooled_within_tissue_metrics.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

pooled_perm <- permutation_test_within_tissue(pooled$actual, pooled$predicted_reference_HRDsum,
                                              pooled$cancer_type, n_perm=1000L)
write.table(pooled_perm, file.path(out_dir,"pooled_within_tissue_permutation.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

cat(sprintf("\nPOOLED  MAE_model=%.3f  MAE_tissue_mean_null=%.3f  skill_vs_tissue_mean=%.3f\n",
            pooled_null$MAE_model, pooled_null$MAE_tissue_mean_null, pooled_null$skill_vs_tissue_mean))
if (is.finite(pooled_null$skill_vs_tissue_mean) && pooled_null$skill_vs_tissue_mean <= 0) {
  cat("WARNING: the model does NOT beat the tissue-mean null. Consistent with the model\n",
      "having learned tissue lineage rather than HRD biology; do not present this as an HRD result.\n")
}

# --- Purity confounding -----------------------------------------------------
if ("purity" %in% names(meta)) {
  pur <- meta$purity[match(pooled$patient_id, meta$patient_id)]
  if (sum(is.finite(pur)) >= 30L) {
    legs <- purity_confound_legs(pooled$actual, pooled$predicted_reference_HRDsum, pur, pooled$cancer_type)
    write.table(legs, file.path(out_dir,"pooled_purity_confound_legs.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    ps <- purity_stratified_metrics(pooled$actual, pooled$predicted_reference_HRDsum, pur,
                                    pooled$cancer_type, mean(meta$HRDsum[!cns]))
    if (!is.null(ps)) write.table(ps, file.path(out_dir,"pooled_purity_stratified.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    pm <- purity_matched_subset(pooled$actual, pooled$predicted_reference_HRDsum, pur,
                                pooled$cancer_type, training_mean=mean(meta$HRDsum[!cns]))
    if (!is.null(pm)) write.table(pm, file.path(out_dir,"pooled_purity_matched.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    cat(sprintf("PURITY  cor(pred,purity)_within=%.3f  cor(label,purity)_within=%.3f\n",
                legs$cor_pred_purity_within, legs$cor_label_purity_within))
  } else {
    cat("Too few finite purity values for stratified analysis; skipping.\n")
  }
} else {
  cat("No purity column in master table; purity confounding NOT assessed.\n")
}

cat("\nmerge complete ->", out_dir, "\n")
