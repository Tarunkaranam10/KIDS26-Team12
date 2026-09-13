# Dataset acquisition and master table

Start with [DATA_ACQUISITION_NOW.md](DATA_ACQUISITION_NOW.md). Labels and six current-GDC beta files were downloaded and checksum verified during setup; this is an engineering subset, not a training cohort. Data are under ignored `data/raw/`; derived tables are under `data/processed/`.

```sh
python scripts/build_master.py
python scripts/prepare_beta.py --smoke-probes 1000
python -m unittest discover -s tests -p "test_*.py"
```

For modeling, replace smoke selection with `--probes config/shared_autosomal_probes.txt`. The deterministic candidate list now contains 384,640 shared autosomal probes with matching hg19 coordinates and `MASK_general=FALSE` in pinned 2022-09 HM450 and EPIC-v1 manifests; see [the bridge specification](17_450K_EPIC_FEATURE_BRIDGE.md). Confirm that PBTP is EPIC v1 before freezing it. For full historical data, run `python scripts/index_publication.py` to derive metadata from its header, followed by `python scripts/build_master.py --metadata config/published_beta_metadata.tsv`. The indexer normalizes mixed patient/sample subtype IDs to patient keys while retaining exact methylation aliquots. Its real 9,664-sample header was checked through a bounded HTTP range request. The team reports that the full matrix has been downloaded, but the 2026-09-13 workspace audit found only the 280,257-byte header fixture under KIDS26; record the actual payload path and checksum before extraction. `--published-matrix` streams beta rows and selects unique specimen columns. Keep engineering and historical outputs in separate directories.

## Matching contract

Patient=12-character TCGA barcode; sample-type key=15 characters; specimen with vial=16; full aliquot retained where supplied. Include primary tumor code 01 only in the initial solid-tumor study. No normal, recurrent, metastatic or unspecified samples. Restrict one methylation file/specimen per patient initially; duplicates and missing metadata go into `matching_audit.tsv`. Do not select replicates by outcome correlation. Resolve later by best predeclared assay QC, a deterministic technical tie-break, and record the adjudication. Never average tumors from different timepoints.

Published HRD has 15-character sample-type IDs. A join to a unique 16-character methylation specimen is therefore a coarser, explicitly flagged match, not proof of the same aliquot. Multiple vials/specimens require manual biospecimen adjudication or exclusion. Separate assays may share patient but not tissue/timepoint. For paired CNV validation, require exact matched specimen where possible and quantify the looser-match sensitivity analysis separately.

The master builder verifies unique HRD rows, nonnegative finite components, sum equality, tumor-type agreement, primary status and published methylation quality exclusions. Missing purity/ploidy stays missing. It excludes when quality annotations are unavailable; do not silently turn unknown QC into pass. Published whitelist QC does not replace raw-array QC.

Required master fields: original file UUID/name/MD5, patient/specimen/aliquot, diagnosis/cancer type, platform/version, preprocessing workflow, LOH/LST/TAI/sum, label source/parameters/build, purity/ploidy and source, matching resolution, QC flags and development/locked-CNS partition. The provided master carries available fields; upstream scar caller/build must be curated before a final scientific freeze.

## Tier gates and storage

Tier 1 labels ~10 MB; full historical betas 41.5 GB; current GDC subset size is printed from live metadata. Tier 2 SNP6+ABSOLUTE segments ~419 MB plus selected raw intensities and normals. Tier 3 RNA ~1.88 GB, optional sequencing much larger and not required today. Record actual bytes, download timestamps, source manifests, software session information and hashes. A checksum verifies bytes, not biological correctness.

| Dataset | Purpose | Expected size | Downloaded? | Verified? | Needed before hackathon? |
|---|---|---:|---|---|---|
| PanCanAtlas/PanImmune labels, subtype, QC, purity/ploidy | Reference HRD labels and joins | 10,249,654 B | Yes | Size, MD5, SHA-256; HRD sum audit | Yes — complete |
| Six current-GDC BRCA beta files | Parser/join engineering | 78,530,106 B | Yes | GDC checksums | No — engineering fixture only |
| Historical PanCanAtlas 450K beta matrix | Main multi-cancer training matrix | 41,541,692,788 B | Team reports complete outside the visible workspace; only header fixture found locally | Header HTTP range verified; payload path/checksum pending | Yes for the planned historical cohort |
| Current-GDC balanced development slice | Small end-to-end current-pipeline engineering bridge | Query-dependent | No | No | Optional; use `--limit-per-project`, never mix silently with historical preprocessing |
| HM450 and EPIC-v1 archived annotation manifests | Deterministic feature bridge | 77,899,903 B | Yes | SHA-256 recorded | Yes — complete, pending PBTP platform confirmation |
| SNP6 plus ABSOLUTE CN resources | Allelic/CN reference and CN benchmark | 419,299,555 B | No | No | Only for gated CN branch |
| Selected masked intensities/IDAT-like GDC inputs | conumee technical pilot | 48,571,534 B for current 3-patient manifest | Manifest only | Manifest query checked | Only for gated CN branch |
| PanCanAtlas RNA | Optional orthogonal comparison | about 1.88 GB | No | No | No |

The 41.5-GB historical matrix is necessary if the team wants the current 7,707-specimen historical cohort. It is not necessary for software smoke testing. `acquire_tcga.py discover --kind beta ... --limit-per-project N` can select a deterministic, balanced current-GDC development slice, but that slice uses a different processing provenance and cannot be pooled with the historical matrix without a technical bridge.

Resume by rerunning the downloader. It refuses a corrupted final file and uses `.part` for incomplete files. GDC-client fallback is documented in the acquisition-now guide. Keep raw download and transformed datasets distinct. Do not mix historical and current-GDC processing until a paired technical bridge is established.

The resolved historical-header audit now retains 7,707/9,664 matrix columns across 32 cancer types. All 8,868 unique matrix participants were resolved to GDC project IDs in `config/published_case_projects.tsv`. The prior 4,397 result was invalidated by incomplete subtype coverage and by a rule that let normal/recurrent files disqualify an otherwise unique primary tumor. The corrected flow counts duplicate candidates among primary tumors only, accepts authoritative GDC project metadata when the subtype table lacks a row, preserves cancer conflicts as exclusions, and reports 40 ambiguous multiple-primary cases. This is a metadata-valid candidate cohort, not a completed beta-level QC cohort. See [cohort QC](16_COHORT_QC.md).
