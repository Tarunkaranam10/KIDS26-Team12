# =============================================================================
# R/provenance.R - Matrix provenance gating (blocker B7).
# =============================================================================
#
# scripts/prepare_beta.py writes a sidecar <matrix_basename>.provenance.json
# next to every beta matrix it builds, recording whether the matrix came from a
# real technical probe allowlist or from the tiny --smoke-probes engineering
# subset. scripts/train_baseline.R has always consulted it; inference did not,
# which meant a fixture could be scored and handed to app/app.R via
# KIDS26_DEMO_RESULTS where it rendered as approved results.
#
# The logic lives here rather than inline so it can be unit tested without a
# fitted model bundle, and so train/inference/merge can converge on one
# convention instead of three near-copies.
#
# DESIGN NOTE - stop vs. warn
# ---------------------------
# Two different failure modes are deliberately treated differently:
#
#   * Fixture or unverifiable matrix -> HARD STOP. There is no legitimate
#     scientific reason to score one, so the caller must say --allow-fixture out
#     loud, and the output is then stamped so the admission travels with it.
#   * Probe allowlist differing from the model's -> WARN ONLY. This is the
#     expected state for genuine transfer work (pediatric EPIC arrays will not
#     share 450k's allowlist), and apply_preprocess() aligns features by name
#     and imputes absentees. Refusing would block the actual research question;
#     silence would hide a platform swap. So: loud, recorded, non-fatal.
# =============================================================================

# Path convention: data/processed/beta.tsv -> data/processed/beta.provenance.json
# Matches the substitution in scripts/train_baseline.R exactly; if one changes,
# change both.
provenance_path <- function(matrix_path) {
  sub("\\.[^.]+$", ".provenance.json", matrix_path)
}

read_matrix_provenance <- function(matrix_path) {
  p <- provenance_path(matrix_path)
  if (!file.exists(p)) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite")
  jsonlite::fromJSON(p)
}

# Classify a sidecar. Returns the provenance class plus a human-readable reason.
# A NULL sidecar is NOT treated as benign: an unlabelled matrix is exactly the
# case this gate exists to catch, since anyone can drop an arbitrary TSV on disk.
classify_matrix_provenance <- function(pr) {
  if (is.null(pr)) {
    return(list(class = "UNVERIFIED_NO_PROVENANCE",
                reason = "no .provenance.json sidecar accompanies this matrix",
                fatal = TRUE))
  }
  if (isTRUE(pr$engineering_only)) {
    return(list(class = "ENGINEERING_FIXTURE_NOT_RESULTS",
                reason = "sidecar declares engineering_only=true (smoke fixture)",
                fatal = TRUE))
  }
  if (is.null(pr$probe_allowlist_sha256)) {
    return(list(class = "UNVERIFIED_NO_ALLOWLIST_HASH",
                reason = "sidecar has no probe_allowlist_sha256 to verify against",
                fatal = TRUE))
  }
  list(class = "verified_research_matrix",
       reason = "sidecar present, non-fixture, allowlist hash recorded",
       fatal = FALSE)
}

# Enforce the classification. `allow_fixture` downgrades a fatal class to a
# recorded admission - it never makes the class disappear.
gate_matrix_provenance <- function(pr, allow_fixture = FALSE) {
  cl <- classify_matrix_provenance(pr)
  if (cl$fatal && !allow_fixture) {
    stop(sprintf(paste0("Provenance gate refused this matrix: %s.\n",
                        "  class: %s\n",
                        "  Scoring it would produce a table whose schema is ",
                        "indistinguishable from real results.\n",
                        "  Pass --allow-fixture to override for engineering ",
                        "work; the output will be stamped as non-scientific."),
                 cl$reason, cl$class), call. = FALSE)
  }
  if (cl$fatal) {
    warning(sprintf("--allow-fixture given: scoring a %s (%s). Output is NOT a scientific result.",
                    cl$class, cl$reason), call. = FALSE)
    cl$class <- paste0(cl$class, ";OVERRIDDEN_BY_ALLOW_FIXTURE")
  }
  cl
}

# Compare the scoring matrix's allowlist against the one the model was trained
# on. train_baseline.R stores the training sidecar as bundle$probe_provenance.
compare_probe_allowlist <- function(model_pr, matrix_pr) {
  a <- if (is.null(model_pr)) NULL else model_pr$probe_allowlist_sha256
  b <- if (is.null(matrix_pr)) NULL else matrix_pr$probe_allowlist_sha256
  if (is.null(a) || is.null(b)) return("unknown_allowlist_not_recorded")
  if (identical(as.character(a), as.character(b))) return("matched")
  "MISMATCH_different_probe_allowlist"
}
