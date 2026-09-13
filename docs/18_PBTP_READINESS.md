# PBTP readiness specification

Protected PBTP inputs stay in approved internal storage. The repository may contain only schemas, code, public/synthetic demonstration fixtures, and approved aggregate outputs. A random project identifier is not itself permission to publish data.

## REQUIRED

| Field or artifact | Specification |
|---|---|
| `project_sample_id` | Random project ID; unique per assay specimen; internal crosswalk stored separately |
| `project_lineage_id` | Random patient/donor lineage shared by patient tumor, PDOX, cell line, and xenograft derivatives |
| `model_type` | Controlled value: `patient_tumor`, `PDOX`, `cell_line`, or `xenograft` |
| tumor/model role | Primary endpoint eligibility: patient tumor versus derived-model secondary analysis |
| methylation assay | Exact 450K, EPIC v1, or EPIC v2 name and manifest version |
| methylation file | IDAT red/green pair preferred; otherwise beta matrix plus exact preprocessing, build, and QC |
| diagnosis | Minimal permitted tumor diagnosis/subgroup sufficient for eligibility and interpretation |
| genomic specimen relationship | Same specimen preferred; otherwise coded same-tumor/timepoint relationship and allowed matching resolution |
| reference `HRDsum` | Independently measured allele-aware genomic value with source/caller/version/build/parameters; required for quantitative validation |
| reference components | HRD-LOH, LST, and TAI values if the supplied HRDsum claims their sum; missing components must remain missing |
| genomic assay QC | WGS/WES/SNP-array type, tumor coverage or assay QC, matched-normal availability, and pass/ambiguity status |
| sharing/use scope | Internal-only/aggregate/demo authorization and responsible data owner |

If an independently measured reference HRDsum cannot be produced, PBTP may receive blinded research predictions, but it cannot serve as labeled quantitative validation.

## STRONGLY DESIRED

| Field or artifact | Use |
|---|---|
| purity and source | Interpret scar calling and methylation/CN signal; missing is unknown |
| ploidy and source | Interpret CN/allelic states and caller behavior |
| allele-specific segments | Chromosome, start, end, total/major/minor copy number, build, coordinate convention |
| methylation intensity QC | Detection, bead/count, sex/concordance, intensity, bisulfite/control, and contamination flags appropriate to the array |
| genomic CNV segments | Total and allele-specific reference for the secondary CN benchmark |
| collection/model timepoint code | Detect temporal or derivative mismatch without exposing direct dates |
| passage number | For PDOX/cell-line/xenograft interpretation, retained internally |
| technical covariates | Batch, plate, center, preservation, extraction, array position where permissible |

## OPTIONAL

| Field or artifact | Use |
|---|---|
| HRDetect/CHORD or other orthogonal score | Contextual concordance; retain its own target definition |
| curated BRCA1/2/RAD51-pathway state | Biological interpretation including biallelic status and reversions |
| BRCA1/RAD51C promoter methylation | Mechanistic annotation; not a substitute for genomic HRDsum |
| treatment/response or functional assay | Exploratory clinical/functional context only under separate approval and endpoint plan |
| RNA or proteomics | Post-hackathon multimodal exploration |

## Grouping and endpoint rules

- The primary external endpoint includes one predeclared eligible **patient tumor per lineage**. Never split samples sharing `project_lineage_id` across training, calibration, or evaluation.
- Derived models are a clustered secondary analysis. Report independent donor count, sample count, within-lineage concordance, and lineage-bootstrap uncertainty.
- Choose among repeated tumors/timepoints using assay and specimen QC defined before outcomes; preserve all exclusions.
- Freeze the feature mask, preprocessing, model, interval/OOD rules, and file hashes before outcome access.
- Do not use PBTP to tune architecture, hyperparameters, feature filters, thresholds, calibration, or abstention. Any adaptation is a separately labeled development analysis requiring another untouched test set.

## Handoff acceptance

The internal custodian should deliver a schema-valid manifest, file checksums, a missingness report, lineage counts, eligible patient-tumor count, derived-model counts, reference-label provenance, and written use scope. The analysis lead returns only approved aggregate metrics and de-identified/precomputed demo outputs to this checkout.
