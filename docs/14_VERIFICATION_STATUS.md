# Verification and handoff status - updated 2026-09-13

## 2026-09-13 readiness update

The team reports that data acquisition is complete. A recursive audit of the shared KIDS26 workspace found the verified label/development/annotation resources described below, but found no file larger than 1 GB and only the 280,257-byte historical-matrix header fixture. Before model-matrix extraction, record the actual 41,541,692,788-byte payload path and verify its checksum. This is a location/provenance gate, not a request to download a second copy.

## Executed successfully

- Python 3.12.14: all Python source files parsed; 17 unit tests passed. Tests cover specimen identity, primary-only duplicate handling, metadata fallback/conflicts, component order/sum, missing/excluded QC, balanced per-project acquisition selection, download corruption/path rejection, CN interval overlap, length-weighted features, centromere boundaries and total-CN non-identifiability of copy-neutral LOH.
- Four PanCanAtlas/PanImmune label/subtype/QC/purity files downloaded: 10,249,654 bytes. All sizes/MD5 verified; SHA256 provenance sidecars retained.
- Six current-GDC open 450K beta files downloaded: 78,530,106 bytes; checksums verified. They are six BRCA patients, not a representative cohort. Matching produced six eligible unique specimens and a 1,000-probe engineering matrix. Current beta format was observed to be headerless two-column data.
- Complete HRD table audited: 10,647 unique sample IDs; zero LOH+LST+TAI sum discrepancies.
- Actual interrupted download resumed from a 1,024-byte partial file and verified. A second run reported verified-existing. The downloader accommodates GDC's observed Content-Range header without the usual bytes prefix while checking the range bounds.
- Historical 450K matrix: bounded 1 MiB HTTP 206 read verified its real header, containing 9,664 sample columns. Only the header was retained as an ignored fixture. Indexing exposed and corrected mixed patient/sample subtype identifiers.
- Historical metadata resolution: all 8,868 unique matrix participants were batch-resolved through the GDC cases API to project IDs. Re-audit retains 7,707/9,664 columns across 32 cancer types, including 642 GBM/LGG locked-CNS candidates. Exclusions in precedence order: 1,363 non-primary, 319 missing/ambiguous HRD, 235 published quality exclusions, and 40 ambiguous multiple-primary candidates. Missing cancer type is now zero; authoritative GDC project metadata fills subtype-table gaps while true conflicts remain excluded. The prior 4,397/22-cancer count is invalid: it was driven by incomplete subtype coverage and a duplicate rule that allowed normal/recurrent files to disqualify a unique primary tumor. The new count is still metadata-level, before full beta extraction and array QC. See `docs/16_COHORT_QC.md`.
- Simulator ran: 17 truth segments and 2,160 noisy probe observations. Its simplified truth is an engineering fixture, not canonical scarHRD validation.
- Literature artifacts contain 70 entries, with 60 mapped local PDFs: 35 previously reviewed, 24 newly available entries from the original list, and one additional chordoma paper absent from the 69-row source. The 25-paper differential review found no exact duplicate or supplementary-only PDF and records explicit absent methods. See `docs/15_NEW_LITERATURE_RECONCILIATION.md`.
- HM450 and EPIC-v1 archived 2022-09 annotation manifests downloaded and SHA-256 verified (77,899,903 bytes combined). Deterministic intersection retains 384,640 shared autosomal `MASK_general=FALSE` CpGs with identical hg19 coordinates; allowlist SHA-256 is `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc`. PBTP EPIC generation remains to be confirmed.
- Public tier manifests generated without downloading the large tiers: historical beta 41,541,692,788 bytes; two CNV files 419,299,555 bytes combined. Current masked-intensity discovery identified 1,586 eligible files; a three-patient/six-file manifest totals 48,571,534 bytes. Intensities were not downloaded. A zero-result raw-intensity query is not proof of universal unavailability.

## Implemented but not executed here

R elastic-net nested LOCO fitting, R smoke tests, frozen inference and Shiny are implemented. Static inspection confirms outer-cancer isolation and training-only variance filtering, imputation, scaling, feature selection, alpha/lambda tuning, null mean and calibration reservation. The training script now exports final coefficients.

**UPDATED 2026-09-17 — the R runtime caveat below is resolved.** The previous
text read: *"R is not installed on Windows or the available WSL image ... so no R
runtime fit, package installation, lockfile, or rendered Shiny validation is
claimed."* R 4.5.0 is available on the cluster, `tests/smoke_model.R` passes
there (jobs 323078428, 323078739), and a full 30-fold nested LOCO fit has now
executed end-to-end (array `323078995`, 30/30 DONE, merged 2026-09-17 06:20).
Leakage assertions pass executably, not just by inspection.

What is now claimed on runtime evidence: nested LOCO fitting, training-only
preprocessing, lambda-path derivation (0/30 folds at a boundary), conformal
interval construction, and the inference provenance gate
(`tests/test_provenance_gate.R`). Results: `docs/24_RESULTS_LOCO_RUN01.md`.

What is still **not** claimed: rendered Shiny validation, a lockfile, any CNS or
pediatric result. The locked CNS partition remains unopened.

The full historical matrix has not been located or streamed end-to-end by this repository. Probe-mask construction is complete for HM450 to EPIC v1, but exact PBTP platform confirmation and empirical preprocessing harmonization remain. Raw-IDAT QC/conumee2 execution, pinned canonical scarHRD comparison, broader platform/allele simulation and PBTP integration remain data-dependent. No biological model has been trained; no transfer accuracy or calibrated HRD-high probability exists.

## First reproducibility commands

Run inside KIDS26-Team12 with Python >=3.10 and R installed as needed:

```sh
python -m unittest discover -s tests -p "test_*.py"
python scripts/acquire_tcga.py published --tier beta --download
python scripts/resolve_tcga_cases.py
python scripts/index_publication.py --case-project-map config/published_case_projects.tsv
python scripts/build_master.py --metadata config/published_beta_metadata.tsv --out data/processed/historical
Rscript scripts/setup.R
Rscript tests/smoke_model.R
```

Inspect the historical matching audit and finalize the technical shared-probe allowlist before extracting model inputs. See the acquisition guide for the complete commands. Keep current-GDC and historical outputs separate.

## Scope and repository hygiene

Original parent plans/papers and team contributor/license files were retained. The template README is preserved verbatim as a historical artifact; its original relative links may no longer resolve from its archive location. The original docs/ai-guidance.md already refers to a nonexistent clone-level AGENTS.md; supplied project instructions reside in the parent directory. Newly authored handoff links were checked. Data, extracted text, models and results are ignored. No Git commit, push or protected-data publication occurred.
