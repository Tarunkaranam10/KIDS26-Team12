# 21. Blocker Resolution Plan

Status date: 2026-09-16. This document consolidates every known blocker into an
ordered, owner-assignable plan. Blockers are grouped by whether they gate the
*first real result* (P0), gate *credible claims* (P1), or gate *deployment and
handoff* (P2). Nothing here asserts a result that has been produced.

The single most important fact: **no biological model has been fit yet.** Every
performance statement in this repository is currently a plan, not a finding.

---

## P0 — Gates the first real result

### B1. Historical matrix payload is unlocated — **RESOLVED 2026-09-16**

**Symptom.** `docs/14_VERIFICATION_STATUS.md` records that only a 280,257-byte
header fixture is present locally, while the expected PanCanAtlas payload is
41,541,692,788 bytes. `config/published_beta_metadata.tsv` carries the expected
md5 (`a92f50490cf4eca98b0d19e10927de9d`).

**Why it blocks.** No feature matrix can be extracted, so no model can be fit.

**Resolution.**
1. Locate the payload on shared scratch; record its absolute path in
   `docs/03_DATA_ACQUISITION.md`.
2. Verify with `md5sum` against the manifest value. A mismatch means re-download,
   not "proceed with caution."
3. Record path + checksum + verification date in `CHANGELOG.md`.
4. Only then run `scripts/prepare_beta.py`.

**Done when.** A checksum-verified path is recorded and `prepare_beta.py` emits a
`.provenance.json` marked as real (not fixture).

> ### RESOLVED 2026-09-16
>
> All four conditions are met. Verified state:
>
> | Item | Value |
> |---|---|
> | Raw payload | `data/raw/pancanatlas/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv` |
> | Size on disk | 41,541,692,788 bytes — **exact match** to the manifest |
> | Expected md5 | `a92f50490cf4eca98b0d19e10927de9d` (`config/published_beta_metadata.tsv`) |
> | sha256 sidecar | `1212f48e8f090fd6afad747e625adf7e66c10a337788aeb81a8a8eaac8f962c4` |
> | Extracted matrix | `data/processed/beta.tsv`, 336,480 probes x 7,707 samples |
> | Provenance | `engineering_only: false`, `track: historical-publication` |
> | Allowlist hash | `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc` |
>
> **How the md5 was verified.** `scripts/acquire_tcga.py` refuses to finalise a
> download unless size AND md5 match the manifest (line 46), and writes the
> `.sha256` sidecar only after that check passes (line 50). The sidecar's
> existence is therefore evidence the md5 gate was satisfied at download time,
> not an independent claim. To re-verify from scratch:
>
> ```bash
> cd data/raw/pancanatlas && sha256sum -c \
>   jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.sha256
> ```
>
> **Note on width.** The extracted matrix has 336,480 probes, not the 384,640 of
> the HM450/EPIC bridge. The bridge is the *allowlist*; 336,480 is its
> intersection with probes actually present in the published matrix. Quote
> 336,480 as the modelling feature count and 384,640 as the platform bridge.
>
> **Still open, and distinct from B1.** The matrix is verified, but per-specimen
> array QC (detection p-values, sex concordance, duplicate audit) has NOT gated
> this cohort — every row carries the same blanket `quality_annotation`. The QC
> gate in `train_baseline.R` now warns about exactly this (see B5).

---

### B2. R runtime never executed — **RESOLVED 2026-09-16**

**Symptom.** `docs/19_R_AND_SHINY_AUDIT.md` states no R executable was available;
all leakage and preprocessing guarantees are verified *by code inspection only*.

**Why it blocks.** Static inspection cannot catch runtime errors, and the
leakage guarantees are the scientific foundation of the whole design.

**Resolution.**
1. Run `Rscript scripts/setup.R` on the cluster; commit `renv.lock` only after a
   successful restore.
2. Run `Rscript tests/smoke_model.R` and record pass/fail in
   `docs/14_VERIFICATION_STATUS.md`.
3. Extend the smoke test with an explicit **leakage assertion**: fit
   preprocessing on a train split, confirm that held-out medians, variance
   ranking, and centering constants are byte-identical to the frozen training
   values and are *not* recomputed on test rows.

**Done when.** Smoke tests pass on real infrastructure and the leakage assertion
is part of the test suite rather than a comment.

> ### SINGLE-FOLD DRY RUN 2026-09-16 (LSF job 323075937, noderome117)
>
> The full modelling path executed end to end on real data for the first time.
> Fold 1 (ACC, the smallest development type) completed successfully:
>
> | Quantity | Value |
> |---|---|
> | Train / test | 6,988 / 77 (locked CNS correctly excluded) |
> | `fit_en` | 53.2 min |
> | Wall time | 55.4 min |
> | **Peak memory** | **173 GB** |
> | Selected hyperparameters | alpha 0.1, lambda 1.411, **interior** to its path |
> | Non-zero coefficients | 772 of 5,000 features |
> | Artefacts written | all five per-fold files |
>
> **Two findings that changed the array configuration.**
>
> 1. **Peak memory is 173 GB, not the 75 GB** the preprocessing-only benchmark
>    implied. `fit_en` additionally holds standardised inner-fold matrices and
>    glmnet's working copies on top of the 19 GB beta matrix. The array was
>    revised from `mem=110GB %15` to `mem=230GB %6`. The original setting would
>    have needed ~3.5 TB against 2.4 TB available and caused mass memory kills.
> 2. **A fold costs ~55 min, not ~45 min**, because the benchmark excluded
>    glmnet. Thirty folds at 6-way concurrency is five waves, roughly 5–6 h.
>
> **The lambda anchoring works as intended.** The selected lambda was interior
> to its data-derived path, which is the outcome that indicates the path was
> wide enough. Under the old hardcoded grid the nearest values were 1 and 10.
>
> **Still outstanding for B2:** `tests/smoke_model.R` has not been run on the
> cluster, and the explicit leakage assertion is not yet in the test suite. The
> dry run exercises the same code path but is not a substitute for that test.

> ### SMOKE TESTS + LEAKAGE ASSERTIONS 2026-09-16 (LSF job 323078428, noderome117)
>
> **B2 is now fully closed.** `tests/smoke_model.R` passed on cluster
> infrastructure, and the leakage guarantee is executable rather than a comment.
>
> Check 7 was added to the suite with six assertions:
>
> | Assertion | What it guards |
> |---|---|
> | 7a | Stored centre/scale equal training-only statistics |
> | 7b | Refitting on the same training rows is deterministic |
> | 7c | Held-out data standardised with training constants stays off-centre |
> | 7d | The transform object is not mutated by being applied |
> | 7e | Features align by NAME with unseen probes present and order shuffled |
> | 7f | The lambda path is reproducible from training rows alone |
>
> **7c is the load-bearing one, and it was validated against a negative
> control.** A deliberately leaky implementation that recentres on the data it
> is handed produces column means of ~7e-17, while the correct implementation
> produces ~1.2. The assertion threshold (>0.5) separates them decisively, so
> this is a test that can actually fail — unlike the tautological QC gate
> described in B5.
>
> One adjustment during development: 7a uses `all.equal(tolerance=1e-12)` rather
> than `identical()`, because `fit_preprocess()` uses `matrixStats::colMeans2`
> while the check uses base `colMeans`. They agree to ~1e-16 — a summation-order
> difference, not a difference in what was computed. Using `identical()` would
> have made the test fail for a reason unrelated to leakage.

---

### B3. Preprocessing is ~42 hours for full nested LOCO — **MITIGATED 2026-09-16**

**Symptom.** `R/model.R` lines 46–62: `fit_preprocess()` uses `apply()` across
~336k columns, ~2.8 minutes per call; nested LOCO needs ~900 calls.

**Resolution.**
1. Replace column-wise `apply()` with `matrixStats::colMedians` / `colVars`
   (documented in-file as reducing the estimate to ~17 hours).
2. Restrict the missingness filter and variance ranking to a single pass.
3. Benchmark on one outer fold before launching the full grid.
4. Submit as an LSF array job, one outer fold per task, rather than a single
   long-running job.

**Caution.** Do not reduce cost by screening features once on the full cohort —
that reintroduces exactly the leakage the design prevents. Speed fixes must be
implementation-level only.

> ### MEASURED 2026-09-16 (LSF job 323075601, nodegpu217)
>
> `matrixStats` reductions are in place and verified numerically identical to
> the old `apply()` path on Inf-free input. Benchmarked on the real
> 336,480 x 7,707 matrix rather than extrapolated:
>
> | Quantity | Value |
> |---|---|
> | `fit_preprocess` | 72.8 s median (73.1 / 71.5 / 72.8) |
> | `apply_preprocess` | 8.7 s |
> | Calls per outer fold | 30 x `fit_preprocess`, 58 x `apply_preprocess` |
> | Per outer fold | ~44.8 min |
> | Serial Phase 1 (30 folds) | **~22.4 h** |
> | Peak RSS | 75 GB |
>
> **An earlier extrapolation from 20k/60k-probe test matrices predicted 19.9 s
> per call and 8.7 h serial. It was wrong by 3.7x.** Cost does not scale
> linearly with probe count once the matrix is ~19 GB, because memory bandwidth
> rather than arithmetic becomes the limit. Treat small-matrix extrapolations in
> this project as lower bounds only.
>
> **Mitigation, not elimination.** Item 4 is what actually makes this tractable:
> `scripts/lsf_loco_array.bsub` runs the 30 folds as independent array tasks at
> 15-way concurrency, giving ~1.5–2 h wall time in two waves. A single
> `fit_preprocess` call is still ~73 s and that has not changed.

**Done when.** One outer fold completes in a measured, recorded wall time and the
full run is scheduled within the event window.

---

### B4. Memory ceiling on matrix load

**Symptom.** `scripts/train_baseline.R` lines 33–72: ~21 GB for `fread()`, 60+ GB
peak through the transpose chain.

**Resolution.** Request a large-memory queue; subset to the 384,640-probe
allowlist *during* read rather than after; avoid retaining intermediate copies;
consider `data.table::setDT` in-place transforms. Record peak RSS.

---

## P1 — Gates credible claims

### B5. QC gate is tautological — **RESOLVED 2026-09-16**

**Symptom.** `scripts/train_baseline.R` line 128 checks
`quality_annotation == "published_450K_no_exclusion"`, but `build_master.py`
hardcodes that exact value. The check can never fail.

**Resolution.** Replace with a real gate: assert the provenance sidecar exists,
that its recorded checksum matches the loaded matrix, and that per-sample QC
fields (detection p-value pass rate, sex concordance) are present and within
bounds. A gate that cannot fail is worse than no gate, because it produces false
assurance in the audit trail.

---

### B6. Latent `sample()` bug — **FIXED 2026-09-16**

**Symptom.** `scripts/train_baseline.R` lines 212–217: `sample(ii, k)` where
`ii` has length 1 silently samples from `1:ii` instead of returning `ii`.

**Why it matters.** Currently masked (smallest development group, OV, has n=10)
but will fire on any subset or rarer stratum, producing invalid fold assignment
without an error.

**Resolution.** Guard with `if (length(ii) == 1L) ii else sample(ii, k)`, or use
`resample <- function(x, ...) x[sample.int(length(x), ...)]`. Add a unit test
with a single-element group.

---

### B7. `predict_frozen.R` has no provenance gate — **RESOLVED 2026-09-16**

**Symptom.** Unlike the training script, inference accepts any beta matrix and
will happily score an engineering fixture, then pass it to the Shiny app via
`KIDS26_DEMO_RESULTS` where it renders as "Approved precomputed results."

**Resolution.** Read the `.provenance.json` sidecar; refuse to score matrices
marked as simulated/fixture unless an explicit `--allow-fixture` flag is passed,
and stamp the output file with the provenance class so downstream display cannot
misrepresent it.

> **Done.** The gate now lives in `R/provenance.R` — shared rather than copied,
> so training and inference cannot drift apart — and runs in
> `scripts/predict_frozen.R` *before* the matrix is read, so an unlabelled 28 GB
> file is refused without being loaded.
>
> Three classes are refused, not one: a declared fixture
> (`engineering_only=true`), a matrix with **no sidecar at all**, and a
> non-fixture matrix whose `probe_allowlist_sha256` is absent. The middle case
> matters most — the original gap was reachable by simply dropping an arbitrary
> TSV on disk, which carries no `engineering_only` flag to catch.
>
> `--allow-fixture` permits engineering runs but cannot launder the label: the
> class is suffixed `;OVERRIDDEN_BY_ALLOW_FIXTURE`, a warning is raised,
> `reportable` is forced `FALSE`, estimates and bounds are blanked, and the
> `provenance` column leads with `NOT A SCIENTIFIC RESULT`. Since `app/app.R`
> validates `KIDS26_DEMO_RESULTS` structurally only and prints `provenance`
> verbatim, putting the class *in that column* is what actually closes the path
> to the demo.
>
> Probe-allowlist mismatch between model and scoring matrix is recorded per row
> (`probe_set_match`) and warned about, but is **not** fatal — cross-platform
> transfer to pediatric EPIC arrays is the research goal, and
> `apply_preprocess()` aligns by feature name.
>
> `tests/test_provenance_gate.R` asserts each refusal executably, including
> against the real `data/processed/beta.provenance.json`.

---

### B8. Lambda grid is hardcoded — **RESOLVED 2026-09-16**

**Symptom.** `R/model.R` uses five fixed lambda values spanning four decades, not
anchored to a glmnet-derived `lambda.max`. A boundary optimum is accepted
silently.

**Resolution.** Derive the path from `glmnet`'s own `lambda.max` per inner fold;
warn loudly if the selected lambda sits at either endpoint of the grid.

> **Done.** `lambda_path()` in `R/model.R` computes
> `lambda.max = max|x'y| / (n * alpha)` per alpha from the fold's own training
> rows, log-spaced down to `0.001 * lambda.max`. Verified to reproduce glmnet's
> internal `lambda.max` exactly at alpha 0.1 / 0.5 / 1.
>
> The old fixed grid is confirmed to have been wasteful: at alpha = 1,
> `lambda.max` is ~2.77, so the hardcoded values 10 and 100 were both guaranteed
> intercept-only fits — two of five grid points were dead.
>
> `check_lambda_boundary()` warns when the winner lands on an endpoint, and the
> fitted bundle records `lambda_paths`, `lambda_at_boundary` and
> `lambda_boundary_side` so the condition survives into the audit trail.
> `scripts/loco_merge.R` reports how many folds hit a boundary and warns if a
> majority did.

---

### B9. `HRD_high_probability` is always NA — **RESOLVED 2026-09-16**

**Symptom.** The exploratory threshold (42) lives in
`config/analysis_protocol.json` but is read by no R code.

**Resolution.** Either (a) drop the column until a calibration set justifies it,
or (b) wire it explicitly and label it exploratory in every output. Option (a) is
preferred before any public presentation — an always-NA column invites a reader
to assume a clinical threshold exists.

> **Done — option (a).** The column was removed from `scripts/predict_frozen.R`
> and replaced with a comment stating why, so the absence is deliberate and
> documented rather than an oversight for someone to "fix" by re-adding it.

---

## P2 — Gates deployment and handoff

### B10. Shiny app governance gap

`app/app.R` validates the schema of `KIDS26_DEMO_RESULTS` but enforces no path
allowlist, checksum, or approval token; raw sample identifiers render as-is; the
`results` object is global and shared across sessions.

**Resolution.** Add a path allowlist + checksum check, hash or alias displayed
identifiers, move `results` inside `server()` for per-session isolation, and make
the provenance label derive from the stamped provenance class rather than from
`nzchar(Sys.getenv(...))`.

### B11. PBTP EPIC generation unconfirmed

The exact EPIC generation is unverified, so the 450K/EPIC bridge cannot be
finalized for pediatric transfer. Confirm generation with the data custodian
before freezing; if EPIC v2, the shared-probe intersection must be rebuilt and
its SHA-256 re-recorded.

### B12. Raw IDAT preprocessing bridge unbuilt

Deferred to post-hackathon per `docs/09_POST_HACKATHON_PLAN.md`. Needed only if
raw IDATs enter the pipeline; currently level-3 betas are used.

---

## Suggested execution order

Status as of 2026-09-16: **B1, B8 and B9 are resolved; B3 and B6 are mitigated.**
The remaining gate before a first real result is B2 (R runtime end-to-end on the
cluster), which the single-fold dry run exercises directly.

```mermaid
graph TD
    B1[B1 locate + checksum matrix<br/>RESOLVED] --> B4[B4 memory-safe load<br/>measured 75 GB peak]
    B2[B2 R runtime + smoke tests<br/>dry run in progress] --> B3[B3 matrixStats + LSF array<br/>MITIGATED]
    B4 --> B3
    B3 --> FIT[First real LOCO fit]
    B5[B5 real QC gate<br/>RESOLVED, warns on blanket QC] --> FIT
    B6[B6 sample bug<br/>FIXED] --> FIT
    FIT --> CONF[Confounding controls: tissue + purity<br/>wired into loco_merge.R]
    CONF --> B7[B7 inference provenance]
    B7 --> B10[B10 app governance]
    B11[B11 PBTP EPIC] --> TRANSFER[Pediatric transfer]

    style B1 fill:#d4edda
    style B5 fill:#d4edda
    style B6 fill:#d4edda
    style B3 fill:#fff3cd
```

The critical path is now: dry run one fold (B2) -> launch the 30-task array ->
`loco_merge.R` -> read the pooled tissue-mean null. **That last step decides
whether there is a result worth reporting at all** (docs/22).
