# Download public TCGA inputs today

Run from the cloned **KIDS26-Team12** repository, with Python >=3.10. These scripts use only the Python standard library. They never request credentials.

```sh
python scripts/acquire_tcga.py published --tier labels --download
python scripts/acquire_tcga.py discover --kind beta --projects TCGA-BRCA TCGA-OV --limit 12 --download
```

The first command obtains published HRD components/sum, subtype annotations, analyte quality annotations and ABSOLUTE purity/ploidy (about 10 MB). The second queries current GDC metadata, selects up to 12 primary-tumor patients, saves UUID/checksum/assay metadata, and downloads their open 450K beta files. It is an engineering subset, NOT a representative performance cohort. Run a separate discovery per cancer if balanced development representation is wanted.

## Tier 1: FAST/MVP

| Resource | Purpose and data level | Acquisition |
|---|---|---|
| TCGA.HRD_withSampleID.txt | Published numeric HRD reference; audit actual component headers and sum | labels command above |
| TCGASubtype.20170308.tsv | Published sample/cancer annotations | labels command |
| merged_sample_quality_annotations.tsv | Whitelist and assay quality adjudication | labels command |
| TCGA_mastercalls.abs_tables_JSedit.fixed.txt | Independent purity/ploidy; QC and confounding audit | labels command |
| Current GDC open 450K beta files | Small development set; current pipeline provenance | discover command above |
| jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv | Historical publication-level 450K matrix; main retrospective experiment | command below |

```sh
python scripts/acquire_tcga.py published --tier beta
python scripts/acquire_tcga.py published --tier beta --download
```

**Storage:** the authoritative manifest reports 41,541,692,788 bytes for the beta matrix (~41.5 GB decimal, ~38.7 GiB). Budget 100-150 GB working disk for this track, depending on conversion and copies; this is a planning allowance, not a measured requirement. Dense 485,000 x 10,000 doubles alone take ~38.8 GB RAM, before copies. Stream/filter by a predeclared probe manifest; use HDF5/DelayedArray for the full matrix. Do not read it with unqualified `read.csv()` into a laptop session. The merged 27K+450K file is NOT an equivalent smaller 450K matrix.

## Tier 2: FULL methylation-CNV benchmark

```sh
python scripts/acquire_tcga.py published --tier cnv --download
python scripts/acquire_tcga.py discover --kind masked-idat --projects TCGA-BRCA --limit 12 --download
```

The SNP6 segmentation + ABSOLUTE allele-specific segment files total about 419 MB. They provide references, NEVER methylation-CNV predictors. Add `--download` to the chosen intensity command only after checking the pairing and selected specimens. Preserve red/green pairs.

**Resolved 2026-09-15 — the earlier zero-result IDAT query was a bug, not a data or permissions problem.** GDC files methylation array intensities *only* as `Masked Intensities`; there is no `Raw Intensities` methylation data type (that value belongs to Affymetrix SNP6 and GeneChip expression arrays, which is why the old `--kind idat` query always returned zero). All 61,830 GDC methylation-array files are **open access**, including 1,790 open 450K IDATs for TCGA-BRCA alone. `--kind idat` and `--kind masked-idat` are now synonyms for the same correct query. No institutional-login automation is needed, because none of this is controlled access.

Masking can remove genotype-sensitive information: masked IDAT is not automatically sufficient for the allele experiment. Verify what the GDC masking step actually removes before relying on it for anything allele-aware.

Methylation-CNV requires IDAT or genuinely retained M/U total intensities with reference controls. Beta ratios alone do not preserve total signal. Obtain assay-compatible normals and test conumee 2.0 on a few paired tumors before expanding. Budget from the manifest byte sum plus 2-3x intermediate space; do not download all raw arrays first.

## Tier 3: optional

```sh
python scripts/acquire_tcga.py published --tier rna --download
```

The v2 expression matrix is ~1.88 GB. Its historical batch correction is not a demonstrated future single-sample transformation. Raw WGS/WES and complex generators are post-hackathon tasks unless independent genomic labels are already prepared.

## Reproducibility and resume

Source manifests in `config/` were obtained from the official [PanImmune](https://gdc.cancer.gov/about-data/publications/panimmune) and [Cell-of-Origin](https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin) resource pages on 2026-09-10. Their UUIDs, expected byte sizes and MD5s are preserved. The downloader streams `.part` files, resumes with HTTP Range, restarts safely if Range is ignored, verifies size+MD5 before rename, and writes SHA256 sidecars. Re-run exactly the same command after interruption. Existing valid files are skipped. A checksum failure stops; inspect the indicated file rather than accepting it.

Manual fallback: open the publication page, choose the exact filenames above, or use the generated filtered manifest:

```sh
gdc-client download -m config/tcga_beta.manifest.tsv -d data/raw/gdc-client
```

Do not download an entire publication manifest by accident. Data Transfer Tool output is UUID-nested; the custom downloader's published output is flat. Retain the layout used and pass paths explicitly downstream. If API access fails, run the same Python commands on an internet-connected research machine; no edits to machine-specific paths are needed.

Current GDC query responses are saved under `data/raw/gdc/discovery`, and flattened metadata/manifests under `config/`. `scripts/build_master.py` and `scripts/prepare_beta.py` provide the next steps (see `03_DATA_ACQUISITION.md`). Review sample matching and quality flags before training.

**Preprocessing boundary:** the historical publication matrix and current SeSAMe GDC beta data are separate tracks. Do not silently mix them. Use the historical matrix for retrospective prediction; establish a bridge using paired processing before claiming production EPIC compatibility. Never normalize PBTP jointly with TCGA to improve external-test agreement.
