# Pre-hackathon checklist and input gates

| Due | Owner role | Required action / output | Acceptance | Fallback |
|---|---|---|---|---|
| Today | Data lead | Run public labels and development beta download; inspect joins | Checksums, sum check, unique specimens, audit saved | Run scripts on connected host |
| Before event | Data lead | Acquire multi-cancer eligible beta/HRD cohort; record cancer Ns and exclusions | At least several represented cancers with enough real samples for outer/inner splits | Narrow to explicit adult proof-of-concept; do not claim pediatric transfer |
| Before event | Methods lead | Confirm PanImmune upstream label provenance and compare a paired reference subset | Same target scale documented; discordance preserved | Use published operational target with explicit provenance limitation |
| Before event | Methylation lead | Confirm PBTP EPIC generation; review generated manifest/build bridge | EPIC-v1 use confirmed or correct version rebuilt; frozen 384,640-probe allowlist hash recorded | Processed 450K-only retrospective scope |
| Before event | Compute lead | Install R/renv/glmnet; run R tests; sufficient disk/RAM | Session report and committed lockfile after validation | Acquisition/data audit + engineering demo only; no untested-model claims |
| Before event | CNV lead | 20-30 paired raw arrays and SNP6/allelic references plus normals | Both channels and assay-matched reference QC | Cut CNV/fusion from hackathon |
| Before event | PBTP custodian | Prepare permitted specimen/lineage manifest and independent genomic reference | Matching, permission and patient-level N known; no PHI in public repo | Defer PBTP; held-out TCGA CNS becomes external-domain experiment |
| Before event | Statistical lead | Freeze endpoints, partitions, metrics, tuning budget and cut rules | Analysis protocol version recorded before test labels revealed | Reduce experiments, preserve test boundary |
| Day 1 start | Team lead | Assign data/model/CNV/review/demo roles | Each deliverable has an owner | Combine roles; cut optional branches |

## PBTP internal input specification

The authoritative required/strongly-desired/optional handoff schema is [docs/18_PBTP_READINESS.md](18_PBTP_READINESS.md). The compact table below remains an operational checklist.

Store protected data outside the public checkout. Use random project IDs (e.g. UUIDs); keep internal crosswalks in approved storage, not this repository or AI inputs. Public files must not contain institutional sample identifiers, PHI, controlled genomic records or unapproved patient predictions. Gitignore is a guardrail, not permission to publish.

| Field / file | Required? | Specification |
|---|---|---|
| project_sample_id | Yes | Random identifier |
| project_lineage_id | Yes | Common patient donor across tumor/PDOX/cell line/xenograft; random ID |
| model_type | Yes | patient_tumor, PDOX, cell_line, xenograft; document passage/timepoint internally |
| diagnosis/subgroup | Yes where permitted | Minimal approved coded subtype, e.g. H3/IDH context |
| methylation | Yes | Red/green IDAT pair preferred; otherwise processed beta plus exact pipeline and QC |
| platform/manifest/build | Yes | Distinguish 450K, EPIC v1, EPIC v2; probe-ID conventions |
| genomic_specimen_match | Yes | Same specimen preferred; date/timepoint concordance coded without direct dates |
| HRDsum/LOH/LST/TAI | For primary quantitative validation | Independent WGS/WES/SNP-based source, caller/version, build, parameters |
| allele-specific segments | Strongly preferred | Chrom/start/end/total/major/minor, coordinate convention |
| purity/ploidy | Strongly preferred | Source and uncertainty; missing is not normal |
| WGS/WES QC | Yes when used as truth | Coverage, capture/assay, matched normal, caller pass/ambiguity |
| optional orthogonal evidence | Optional | HRDetect/CHORD with QC; curated biallelic defects; functional assays |
| permitted sharing scope | Yes | Internal-only vs approved public aggregate/demo; owner authorization |

Primary PBTP endpoint uses patient tumors, one predeclared specimen per donor. Derived models are a clustered secondary analysis with within-lineage comparisons. Count independent donors in uncertainty and all splits. No automatic internal data upload, credential use or Git push is included.

## Go/no-go

Do not start expensive model tuning until the target audit and group identities pass. Do not claim EPIC deployability until the preprocessing bridge is checked. Do not expose PBTP outcomes before the model/configuration hash is recorded. R execution, full-dataset acquisition and PBTP input preparation must be completed by the team before a real-data model run.
