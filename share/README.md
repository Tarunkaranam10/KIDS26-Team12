# KIDS26 Team 12 — shared public TCGA dataset

**Shared by:** Evan Savage (esavage), St. Jude KIDS26 Biohackathon Team 12
**Location:** `/lustre_scratch/shared_scratch/kids26_team12_share`
**Staged:** 2026-09-15
**Source repository:** `KIDS26-Team12` (`/research/groups/shelagrp/home/esavage/KIDS26/KIDS26-Team12`)
**Approximate total size:** ~69 GB

> **Point-in-time staging copy.** This directory and any documentation copied into it
> reflect the repository as of the staging date above. It is not kept in sync with the
> source repository; for current documentation and code, use the `KIDS26-Team12` repository.

---

## 1. What this is

This directory contains the **public, open-access TCGA/PanCanAtlas inputs** assembled for Team 12's
project: *Can DNA methylation predict independently measured homologous-recombination-deficiency
(HRD) genomic-scar burden across cancer types?*

Everything here is **public data** downloaded from the NCI Genomic Data Commons (GDC) and the
PanCanAtlas/PanImmune and Cell-of-Origin publication resource pages. There are **no protected,
controlled-access, or patient-identifiable data** in this share. No PBTP (pediatric brain tumor
program) data, no St. Jude internal data, and no controlled-access GDC files are included.

Alongside the data we include the acquisition/preparation **scripts, manifests, and documentation**,
so you can verify provenance and regenerate every derived file.

### Scientific framing (important for interpretation)

Canonical HRDsum = HRD-LOH + LST + TAI, and it is measured from **allele-aware genomic data**
(SNP6 arrays / sequencing). Methylation arrays do not directly observe loss of heterozygosity or
telomeric allelic imbalance. In this project methylation is therefore treated as a **predictor of an
independently measured continuous reference HRDsum**, not as a way to re-derive the scar score.
Nothing here is a validated clinical HRD assay.

---

## 2. Directory layout

```
kids26_team12_share/
├── README.md                     <- this file
├── repo_README.md                <- technical README from the source repository
├── repo_README_BIOHACKATHON.md   <- project profile, mission, 3-day plan
├── config/                       <- manifests, probe bridge, metadata crosswalks (~15 MB)
├── docs/                         <- full project documentation set
├── scripts/                      <- acquisition and preparation code (Python + R)
└── data/
    ├── raw/
    │   ├── pancanatlas/          <- published PanCanAtlas matrices and label tables (~42 GB)
    │   ├── annotation/           <- HM450 / EPIC-v1 probe manifests (75 MB)
    │   └── gdc/
    │       ├── beta/             <- 13 current-GDC SeSAMe level-3 beta files (166 MB)
    │       ├── masked-idat/      <- 24 open 450K IDATs, 12 BRCA tumors, Red/Grn pairs (187 MB)
    │       └── discovery/        <- raw GDC API query + response JSON (provenance)
    └── processed/                <- derived analysis-ready tables (~27 GB)
```

---

## 3. Raw data inventory (`data/raw/`)

Every file in `data/raw/pancanatlas/` and `data/raw/gdc/beta/` has a `.sha256` sidecar written at
download time, **after** size and MD5 were verified against the official source manifest.

### 3.1 `data/raw/pancanatlas/` — published PanCanAtlas resources

| File | Bytes | Tier | What it is |
|---|---:|---|---|
| `jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv` | 41,541,692,788 | 1 | **Primary predictor matrix.** Historical PanCanAtlas Illumina HumanMethylation450 beta-value matrix, whitelisted aliquots. ~485K probe rows × 9,664 sample columns. Values are beta in [0,1]; missing values present. |
| `TCGA.HRD_withSampleID.txt` | 280,336 | 1 | **Primary label source.** Published HRD components (HRD-LOH, LST, TAI) and HRDsum, keyed by 15-character TCGA sample-type ID. 10,647 unique sample IDs; component sum verified consistent with the reported sum. |
| `TCGASubtype.20170308.tsv` | 603,836 | 1 | Published sample/cancer-type and molecular-subtype annotations. Uses a mix of patient- and sample-level TCGA IDs (normalized by our indexer). |
| `merged_sample_quality_annotations.tsv` | 8,463,670 | 1 | PanCanAtlas whitelist / assay-quality adjudication. Used to exclude aliquots with "Do not use" style annotations. |
| `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | 901,812 | 1 | ABSOLUTE purity and ploidy per sample. Independent QC and confounding covariate. |
| `TCGA_mastercalls.abs_segtabs.fixed.txt` | 253,061,161 | 2 | ABSOLUTE allele-specific copy-number segment table. **Reference/benchmark only** — never a methylation-derived predictor. |
| `broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg` | 166,238,394 | 2 | Affymetrix SNP6 genome-wide segmentation (whitelisted). The allele-aware substrate from which canonical scar scores derive. Reference only. |
| `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv` | 1,879,492,443 | 3 | **Tier 3, optional.** Batch-adjusted PanCanAtlas RNA-seq v2 gene expression matrix. Included for orthogonal comparison. See caveat in §7. |

SHA-256 (also in the `.sha256` sidecars):

```
1212f48e8f090fd6afad747e625adf7e66c10a337788aeb81a8a8eaac8f962c4  jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv
35af0e8811dd47fd9954b052a0bc9900ac128f3fcf83ec6bcbdc0a7fb86e3e8d  TCGA.HRD_withSampleID.txt
1631cf7dd1324074957b42ee2ecd556b24d0977e5b2b2fb10c891c641994eb94  TCGASubtype.20170308.tsv
e3190253ccb55f9f5ccc558b7f994577b2e90c2ee306814231b7d5f077364deb  merged_sample_quality_annotations.tsv
f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd  TCGA_mastercalls.abs_tables_JSedit.fixed.txt
04080d7705471dc6988326cb64457c9e59bad7e9af3d7562ae9a7a615ab238ea  TCGA_mastercalls.abs_segtabs.fixed.txt
dd27bcf0068a537ec052653d4a856a85e5e47b2471196e846dee138758be4504  broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg
d455f4be62a65b282974e74dbffda94d31d35a5c238ecb32c625dc459883f5cf  EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv
```

### 3.2 `data/raw/annotation/` — probe manifests

| File | Bytes | What it is |
|---|---:|---|
| `HM450.hg19.manifest.202209.tsv.gz` | 28,670,538 | Pinned 2022-09 HM450 probe manifest, hg19 coordinates, with `MASK_general` flags |
| `EPIC.hg19.manifest.202209.tsv.gz` | 49,229,365 | Pinned 2022-09 EPIC-v1 probe manifest, hg19 coordinates |

These two define the **450K↔EPIC feature bridge**: probes present on both platforms, with matching
hg19 coordinates and `MASK_general = FALSE`. That intersection is 384,640 probes (see §4).

### 3.3 `data/raw/gdc/` — current-GDC engineering subset

`beta/` holds **13 open-access SeSAMe level-3 beta files** (`*.methylation_array.sesame.level3betas.txt`),
one directory per GDC file UUID, from TCGA-BRCA and TCGA-OV primary tumors. Total ~166 MB.

> **These are an engineering fixture, not a cohort.** They exist so that parsers, joins, and the
> end-to-end pipeline can be exercised on current-GDC-provenance files. They have **different
> preprocessing** (SeSAMe, current GDC workflow) than the historical publication matrix and must not
> be pooled with it without an explicit paired technical bridge.

`discovery/` holds the literal GDC API request and response JSON for each discovery query
(`beta_*`, `idat_*`, `masked-idat_*`), preserved so the query is auditable and reproducible.

**Note on raw intensities (IDAT):** an earlier discovery run returned zero files and was initially
recorded as an open question. **It was a query bug, now fixed (2026-09-15).** GDC files methylation
array intensities *only* under the data type `Masked Intensities`; there is no `Raw Intensities`
methylation data type, so the old query could never match. All 61,830 GDC methylation-array files
are **open access**, including 1,790 open 450K IDATs for TCGA-BRCA.

`masked-idat/` therefore now holds **24 open-access 450K IDAT files** (~187 MB): 12 TCGA-BRCA
primary tumors, each with its **Red and Green pair intact** (12 Red + 12 Grn, verified). File
magic bytes confirm genuine `IDAT` format. This is the input for the conumee 2.0 methylation-CNV
pilot. See §7 item 4 for what is still required before that pilot means anything.

---

## 4. Config and metadata (`config/`)

| File | What it is |
|---|---|
| `panimmune_source_manifest.tsv`, `celloforigin_source_manifest.tsv` | Official source manifests (UUID, filename, MD5, size) retrieved 2026-09-10 from the GDC PanImmune and Cell-of-Origin publication pages. The root of provenance for everything in `data/raw/pancanatlas/`. |
| `tcga_labels.manifest.tsv`, `tcga_beta.manifest.tsv`, `tcga_cnv.manifest.tsv`, `tcga_rna.manifest.tsv` | Filtered per-tier manifests generated by the acquisition script; also usable directly with `gdc-client download -m <manifest>`. |
| `gdc_beta_all.tsv` / `gdc_beta_dev.tsv` / `gdc_beta_dev.manifest.tsv` | Flattened current-GDC beta discovery results: all eligible files, the selected development subset, and a gdc-client manifest. |
| `gdc_idat_*.tsv`, `gdc_masked-idat_*.tsv` | Same, for raw and masked intensities. `gdc_idat_all.tsv` is header-only (zero hits). |
| `shared_autosomal_probes.txt` | **The frozen feature bridge.** 384,640 autosomal probes shared by HM450 and EPIC-v1 with matching hg19 coordinates and `MASK_general=FALSE`. SHA-256 `f359e43bc171c544e55e60000905bb0bba1c44ee0b37f3144c8ec733e052e5fc`. |
| `shared_autosomal_probes.metadata.json` | Build parameters and manifest versions behind the probe list. |
| `published_beta_metadata.tsv` | Sample metadata derived from the historical matrix header (aliquot → patient/sample/cancer type). |
| `published_case_projects.tsv` | GDC-authoritative participant → project-ID crosswalk. All 8,868 unique matrix participants resolved. |
| `analysis_protocol.json` | Predeclared analysis protocol (splits, preprocessing, model family). |
| `literature_*.json`, `*_extras.json` | Curated literature-evidence records backing the method choices. |

---

## 5. Processed data (`data/processed/`)

These are **derived** from the raw files by the scripts in `scripts/`. Regenerate rather than edit.

| File | Size | Description |
|---|---:|---|
| `beta.tsv` | ~27 GB | Analysis-ready beta matrix from the **historical publication track**. Tab-separated, first column `probe_id`, then 7,707 specimen columns. **336,480 probe rows** (the 384,640-probe bridge intersected with probes actually present and retained in the matrix). |
| `beta.provenance.json` | 512 B | Provenance for the above: `n_samples=7707`, `n_probes=336480`, `track="historical-publication"`, probe-allowlist SHA-256, and the matrix SHA-256 `e3642d3038b3298f302c59ce12801efa5ce8aead4fab470756eab7f46290062a`. |
| `master_samples.tsv` | 3.4 MB | **The master sample table** — 7,707 rows, one per eligible specimen. This is the file to join against. Columns below. |
| `matching_audit.tsv` | 3.3 MB | 9,664 rows — every candidate matrix column with its disposition (kept, duplicate, wrong sample type, QC-excluded, metadata-missing). **Read this before trusting the cohort**; exclusions are explicit here, not silent. |

### 5.1 `master_samples.tsv` columns

`id`, `filename`, `md5`, `size`, `patient_id`, `sample_id`, `cancer_type`, `cancer_type_source`,
`platform`, `workflow`, `access`, `HRD_LOH`, `LST`, `TAI`, `HRDsum`, `purity`, `ploidy`,
`label_source`, `match_resolution`, `quality_annotation`, `partition`

### 5.2 Cohort composition (n = 7,707 primary tumors, 32 cancer types)

| Partition | n | Meaning |
|---|---:|---|
| `development` | 7,065 | Available for model development / leave-one-cancer-out CV |
| `locked_CNS` | 642 | GBM + LGG, **held out** as the adult-CNS proxy for pediatric transfer. Do not touch during development. |

Per-cancer counts:

| Cancer | n | Cancer | n | Cancer | n | Cancer | n |
|---|---:|---|---:|---|---:|---|---:|
| BRCA | 743 | HNSC | 509 | LGG | 507 | THCA | 464 |
| UCEC | 403 | LUAD | 394 | BLCA | 388 | PRAD | 387 |
| STAD | 370 | LIHC | 355 | LUSC | 355 | KIRC | 297 |
| CESC | 290 | COAD | 278 | KIRP | 262 | SARC | 213 |
| PCPG | 160 | ESCA | 154 | TGCT | 149 | PAAD | 143 |
| GBM | 135 | THYM | 107 | SKCM | 104 | READ | 92 |
| UVM | 80 | MESO | 78 | ACC | 77 | KICH | 65 |
| UCS | 56 | DLBC | 47 | CHOL | 35 | OV | 10 |

> **OV n=10 is a real limitation.** Ovarian carcinoma is the canonical HRD-enriched cancer type, and
> it is almost absent from this cohort because most TCGA-OV methylation was assayed on the 27K
> platform rather than 450K. Any claim about HRD-high biology should account for this.

### 5.3 Matching contract (how specimens were joined)

- Patient key = 12-character TCGA barcode; sample-type key = 15 characters; specimen with vial = 16.
- **Primary tumor (code `01`) only.** No normal, recurrent, metastatic, or unspecified samples.
- One methylation specimen per patient; duplicates and metadata gaps go to `matching_audit.tsv`.
- All 7,707 rows carry `match_resolution = sample_type_15char_unique_methylation_specimen`.

**Read this carefully:** published HRD labels are keyed at 15 characters, while methylation
specimens are 16 characters. Every join in this cohort is therefore a **coarser sample-type-level
match, explicitly flagged — not proof of the same physical aliquot.** For any analysis where
aliquot identity matters (especially paired CNV validation), treat this as a known source of
mismatch and quantify it.

Additional guarantees enforced by the builder: unique HRD rows, non-negative finite components,
HRD component sum equals reported sum, tumor-type agreement, primary status, and published
methylation-quality exclusions. Missing purity/ploidy stays missing. Specimens with **unavailable**
quality annotations are excluded rather than assumed to pass.

---

## 6. How to reproduce / extend

All acquisition scripts are **Python ≥3.10, standard library only**, and never request credentials.
On this cluster use a conda Python (e.g. `/home/esavage/miniconda3/bin/python`) — the default
`python` on some login nodes is 3.5.2 and will fail with a syntax error on f-strings.

```sh
# Re-download and re-verify any tier (existing verified files are skipped)
python scripts/acquire_tcga.py published --tier labels --download
python scripts/acquire_tcga.py published --tier beta   --download   # 41.5 GB
python scripts/acquire_tcga.py published --tier cnv    --download
python scripts/acquire_tcga.py published --tier rna    --download

# Re-run current-GDC discovery
python scripts/acquire_tcga.py discover --kind beta --projects TCGA-BRCA TCGA-OV --limit 12

# Rebuild derived tables
python scripts/index_publication.py
python scripts/build_master.py --metadata config/published_beta_metadata.tsv
python scripts/prepare_beta.py \
  --published-matrix data/raw/pancanatlas/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv \
  --probes config/shared_autosomal_probes.txt
```

The downloader streams to `.part` files, resumes with HTTP Range, verifies size + MD5 before
renaming, and writes a SHA-256 sidecar. Re-running the same command after an interruption is safe
and idempotent. A checksum failure **stops** rather than accepting the file.

### Verifying this copy yourself

```sh
cd /lustre_scratch/shared_scratch/kids26_team12_share/data/raw/pancanatlas
cat *.sha256 | sha256sum -c -
```

Historical note: four sidecars were briefly CRLF, because `Path.write_text()` translates `\n` to
`os.linesep` and the original download ran on Windows. That affected the sidecar text only, never
the data bytes. Fixed at the source in `scripts/acquire_tcga.py` (LF is now pinned explicitly);
all sidecars in this share are LF and verify cleanly without preprocessing.

(The staging job runs this check itself; output is in the job log.)

---

## 7. Caveats, limits, and things not to do

1. **Memory.** A dense 485,000 × 10,000 double matrix is ~38.8 GB before copies. Do not load
   `beta.tsv` or the raw 41.5 GB matrix with an unqualified `read.csv()` / `pd.read_csv()`.
   Stream by rows, subset by a predeclared probe list, or use HDF5/DelayedArray.
2. **Two preprocessing tracks must not be mixed.** The historical publication matrix
   (`data/raw/pancanatlas/...betaValue_whitelisted.tsv` → `data/processed/beta.tsv`) and the current
   GDC SeSAMe files (`data/raw/gdc/beta/`) have different provenance. Pooling them without a paired
   technical bridge will manufacture batch effects that look like signal.
3. **CNV/segment files are references, never methylation predictors.** `*abs_segtabs*` and the SNP6
   `.seg` are the allele-aware ground-truth substrate. Using them as model inputs is circular.
4. **Beta values cannot yield copy number.** Beta is a ratio and does not preserve total intensity.
   Methylation-CNV requires IDAT or genuinely retained M/U intensities with reference controls.
   The 24 IDATs in `data/raw/gdc/masked-idat/` are a **12-tumor pilot input, not a CNV benchmark**.
   Before they support any claim you still need assay-compatible **normal reference** arrays, and
   you must confirm what the GDC masking step removes — masked IDATs are not automatically
   sufficient for anything allele-aware (which includes the LOH and TAI components of HRDsum).
5. **The RNA matrix is historically batch-corrected.** That correction was fit on the whole cohort;
   it is not a demonstrated single-sample transformation, so it is not directly applicable to a
   prospective new sample.
6. **`locked_CNS` (n=642) is held out by design.** Using it during development destroys the only
   adult-CNS proxy for the pediatric transfer question.
7. **A checksum verifies bytes, not biological correctness.** This is a metadata-valid candidate
   cohort. Beta-level QC (detection p-values, array-level outliers, sex checks) has **not** been
   completed. See `docs/16_COHORT_QC.md`.
8. **Not a clinical assay.** No result derived from these data is a validated HRD test, a functional
   repair assay, or a treatment recommendation.

---

## 8. Provenance and licensing

Source manifests were retrieved 2026-09-10 from the official GDC publication resource pages:

- PanImmune: https://gdc.cancer.gov/about-data/publications/panimmune
- PanCan Cell-of-Origin: https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin

All files are **open-access TCGA data**. Use is governed by the NCI GDC data use policies and the
TCGA publication guidelines; please cite the originating TCGA/PanCanAtlas publications for any
derived work. The Team 12 code in `scripts/` carries the repository license (`LICENSE.md` in the
source repo).

---

## 9. Contact

Questions about the cohort construction, the matching contract, or anything that looks wrong:
**Evan Savage (esavage)**, KIDS26 Team 12. Further documentation is in `docs/`, starting with
`docs/DATA_ACQUISITION_NOW.md` (acquisition tiers) and `docs/03_DATA_ACQUISITION.md`
(cohort and master-table construction).

If something here disagrees with `docs/`, trust `data/processed/beta.provenance.json` and the
`.sha256` sidecars — those were written by the code at the time the bytes landed.
