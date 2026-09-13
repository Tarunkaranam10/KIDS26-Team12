# Decisions and changes

## 2026-09-10 - independent methods review and prototype scaffolding

- Inspected parent KIDS26 and cloned KIDS26-Team12; preserved source papers, prior plans, contributor/license/team files, and original template README.
- Selected continuous reference-HRDsum methylation elastic net; gated CNV-only/fusion comparisons. Deferred deep generative models and RNA until core evidence exists.
- Distinguished allele-specific canonical scars from total-CN proxies and from current repair function/therapy response.
- Corrected interpretation of Chen MLPA target, Rempel SNP6 pan-cancer truth, DiffuCpG inpainting, SyntheVAEiser expression fine-tuning and sample-specific classifier retraining.
- Added 69-row evidence CSV/JSON and human digest; 35 local PDFs reviewed by relevant methods/pages, 34 missing-full-text entries left explicit. Source-year/version discrepancies retained.
- Added staged download, UUID/MD5 manifests, conservative specimen joins, historical-header indexer, streaming/bounded beta preparation and ignored data directories.
- Verified small real public downloads, six matched BRCA samples, sum integrity of all 10,647 published HRD rows, actual HTTP resume, and Python engineering tests.
- Added R nested baseline, frozen inference, targeted R tests and precomputed Shiny demo. No R executable available: R execution and package lock creation remain pending.
- Added formal checklist, runbook, validation/calibration/OOD rules, PBTP governance spec, risks and research roadmap.
- No PBTP data accessed, no institutional-login automation, no commit/push.

## 2026-09-10 - differential literature reconciliation and readiness audit

- Reconciled 25 newly available PDFs against the frozen 35-PDF inventory: 24 fill original map rows and one chordoma paper is additional evidence. Expanded the evidence matrix to 70 entries/60 local PDFs.
- Primary architecture unchanged: shared-probe methylation elastic-net regression to continuous allele-aware reference HRDsum. New TAI/HRDsum/GSS papers strengthen the non-identifiability of canonical LOH/TAI from ordinary total CN.
- Refined implementation using new consensus/classifier evidence: explicit label caller/build/parameters and tumor-specific thresholds, a versioned 450K/EPIC-v1 bridge, per-sample coverage/OOD abstention, and patient-lineage grouping.
- Resolved all historical matrix participants to GDC projects and corrected the provisional cohort from an invalid 4,397/22-cancer result to 7,707 metadata-eligible columns across 32 cancers. Preserved all exclusions in cohort-QC artifacts.
- Added deterministic annotation acquisition/bridge construction, balanced development-download selection, final coefficient export, PBTP readiness specification, and precomputed Shiny validation summaries.
- Synthetic augmentation remains deferred; no new generator study provides external methylation-HRD transfer evidence.
- No biological model trained, PBTP data accessed, commit, or push.

## 2026-09-13 - Biohackathon project profile and team coordination

- Added a completed Biohackathon-facing project profile with the scientific question, inputs, outputs, stack, mission, three-day milestones, limitations, and navigation.
- Filled the project plan, team roles, checklist handoff, Slack channels, mentor contact, meeting notes, and reusable check-in record using the supplied team information.
- Preserved the detailed technical README and the original blank organizer template; the archive now points readers to the completed version.
- Recorded the team's statement that data are downloaded while retaining the reproducibility finding that the 41.5-GB historical matrix is not visible under the KIDS26 workspace and still needs a confirmed path/checksum before extraction.
