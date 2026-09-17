# =============================================================================
# tests/test_provenance_gate.R - Executable assertions for blocker B7.
# =============================================================================
#   Rscript tests/test_provenance_gate.R     (from the repository root)
#
# These are refusal tests. Each one asserts that a specific bad input is
# REJECTED, because the failure mode B7 describes is a gate that exists but
# does not actually stop anything.
# =============================================================================

source("R/provenance.R")

real <- list(engineering_only = FALSE, n_samples = 7707L, n_probes = 336480L,
             probe_allowlist_sha256 = "f359e43bc171c544e55e60000905bb0bba1c44ee")
fixture <- list(engineering_only = TRUE, n_samples = 8L, n_probes = 50L,
                probe_allowlist_sha256 = "deadbeef")
nohash <- list(engineering_only = FALSE, n_samples = 100L, n_probes = 1000L)

# 1. The real matrix passes and is classed as such.
stopifnot(identical(gate_matrix_provenance(real)$class, "verified_research_matrix"))

# 2. An engineering fixture is refused outright.
stopifnot(inherits(try(gate_matrix_provenance(fixture), silent = TRUE), "try-error"))

# 3. A matrix with NO sidecar at all is refused. This is the case an attacker or
#    a hurried teammate hits: drop any TSV on disk with no metadata beside it.
stopifnot(inherits(try(gate_matrix_provenance(NULL), silent = TRUE), "try-error"))

# 4. A non-fixture matrix whose allowlist cannot be verified is still refused -
#    "not declared fake" is not the same as "verified real".
stopifnot(inherits(try(gate_matrix_provenance(nohash), silent = TRUE), "try-error"))

# 5. --allow-fixture permits the run but must NOT launder the label: the class
#    still names the fixture and additionally records that it was overridden.
ov <- suppressWarnings(gate_matrix_provenance(fixture, allow_fixture = TRUE))
stopifnot(grepl("ENGINEERING_FIXTURE", ov$class),
          grepl("OVERRIDDEN_BY_ALLOW_FIXTURE", ov$class))

# 6. The override must also raise a warning, not pass quietly.
stopifnot(length(withCallingHandlers(
  { w <- character(0)
    gate_matrix_provenance(fixture, allow_fixture = TRUE); w },
  warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") }
)) > 0)

# 7. Allowlist comparison distinguishes all three states. A mismatch is a
#    warning-level condition, not a refusal, because cross-platform transfer is
#    the actual research goal (see R/provenance.R design note).
stopifnot(identical(compare_probe_allowlist(real, real), "matched"),
          identical(compare_probe_allowlist(real, fixture), "MISMATCH_different_probe_allowlist"),
          identical(compare_probe_allowlist(real, nohash), "unknown_allowlist_not_recorded"),
          identical(compare_probe_allowlist(NULL, real), "unknown_allowlist_not_recorded"))

# 8. The sidecar path convention must match the one train_baseline.R uses.
stopifnot(identical(provenance_path("data/processed/beta.tsv"),
                    "data/processed/beta.provenance.json"))

# 9. End-to-end against the real sidecar on disk, if present: the production
#    matrix must pass its own gate. Guards against the sidecar drifting out of
#    the schema the gate expects.
if (file.exists("data/processed/beta.provenance.json")) {
  stopifnot(identical(
    gate_matrix_provenance(read_matrix_provenance("data/processed/beta.tsv"))$class,
    "verified_research_matrix"))
  cat("  (verified against the real data/processed/beta.provenance.json)\n")
}

cat("B7 provenance gate tests passed: fixtures, unlabelled matrices and\n")
cat("  unverifiable allowlists are all refused; overrides are stamped.\n")
