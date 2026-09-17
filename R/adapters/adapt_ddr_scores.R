# adapt_ddr_scores.R
# Reshapes beta_hrd_genes.tsv (wide: probes x samples) into long format,
# joins gene names, and joins sample metadata/HRD scores.
# Output: one tidy data.frame ready for VizModules.
#
# PATHS: these files live in front_end_data/ddr_scars/, not data/ddr_scars/.
# data/ is gitignored cohort storage; front_end_data/ is the small committed
# extract prepared for the front-end team.
#
# CAUTION: master_samples.tsv carries a `partition` column that encodes the
# locked GBM/LGG split, plus per-sample purity and ploidy. Those are analysis
# covariates and partition bookkeeping - do not feed them to a model or use
# them to filter what gets displayed as a result.

library(data.table)

adapt_ddr_scores <- function(
  data_dir     = file.path("front_end_data", "ddr_scars"),
  beta_path    = file.path(data_dir, "beta_hrd_genes.tsv"),
  probes_path  = file.path(data_dir, "hrd_gene_probes.tsv"),
  samples_path = file.path(data_dir, "master_samples.tsv")
) {
  # Fail with the missing path rather than a bare fread() error, since the
  # original default pointed at a directory that does not exist.
  paths <- c(beta = beta_path, probes = probes_path, samples = samples_path)
  absent <- paths[!file.exists(paths)]
  if (length(absent)) {
    stop("Missing input file(s):\n  ", paste(absent, collapse = "\n  "),
         "\nRun from the repository root, or pass data_dir explicitly.",
         call. = FALSE)
  }

  beta <- fread(beta_path)
  # hrd_gene_probes.tsv is headerless: probe_id then gene symbol.
  probes <- fread(probes_path, header = FALSE, col.names = c("probe_id", "gene"))
  samples <- fread(samples_path)

  # wide -> long: one row per probe x sample
  beta_long <- melt(
    beta,
    id.vars = "probe_id",
    variable.name = "sample_id",
    value.name = "beta_value"
  )

  # melt() returns sample_id as a factor over the column names; the metadata
  # table stores it as character. Coerce so the join keys match on type.
  beta_long[, sample_id := as.character(sample_id)]

  # add gene name per probe
  beta_long <- merge(beta_long, probes, by = "probe_id", all.x = TRUE)

  # join sample metadata + HRD scores
  out <- merge(beta_long, samples, by = "sample_id", all.x = TRUE)

  # Surface join failures instead of letting them travel downstream as NA.
  unmatched <- sum(is.na(out$HRDsum))
  if (unmatched) {
    warning(unmatched, " of ", nrow(out),
            " probe x sample rows had no metadata match.", call. = FALSE)
  }

  out[]
}

# Example use:
# ddr_data <- adapt_ddr_scores()
# head(ddr_data)
