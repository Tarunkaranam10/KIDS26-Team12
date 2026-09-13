# Additional primary-source checks (2026-09-10)

These supplement the local evidence matrix. A webpage lookup is not a complete supplementary-methods audit.

| Source | Evidence checked | KIDS26 implication / remaining gap |
|---|---|---|
| [scarHRD repository](https://github.com/sztup/scarHRD) and [entry point](https://raw.githubusercontent.com/sztup/scarHRD/master/R/scar_score.R) | Allele-specific input columns; build and LOH cutoff; component calls and unadjusted returned sum | Pin revision and chromosome tables; full package execution pending |
| [Sztupinszki 2018 paired SNP6/WES figure](https://discovery.ucl.ac.uk/10053224/1/Migrating%20the%20SNP%20array-based%20homologous%20recombination%20deficiency%20measures%20to%20next%20generation%20sequencing%20data%20of%20breast%20cancer.pdf) | Component r .73-.84; sum .87 | Correlation supports feasibility, not assay interchangeability |
| [Popova 2012](https://aacrjournals.org/cancerres/article-abstract/72/21/5454/576090) | >=10Mb adjacent regions; <3Mb filtering; arm-wise counts | Centromeres and smoothing matter |
| [Birkbak 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/) | Telomeric AI excludes centromere-crossing regions | Total CN telomere burden is a different quantity |
| [conumee2](https://pmc.ncbi.nlm.nih.gov/articles/PMC10868300/) | TCGA LUSC n367 SNP comparison r .91; LGG focal validation n239 | Total-CN feasibility; raw caller execution not performed here |
| [Zhou et al. 2017, PDF pp.1-2](https://zwdzwd.s3.amazonaws.com/papers/2017NAR.pdf) | 59 explicit +1052 extra SNP-informative EPIC loci; probe masks/channel-switch chemistry; matched WGBS/450K and normals | Genotyping utility is not evidence of tumor LOH/TAI reconstruction; exact modern manifest still required |
| [expHRD 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11241885/) | Elastic net then bootstrap 356-gene signature and ssGSEA; random TCGA test; OV/GDC evaluations | Useful n-of-1 comparator; GDC label does not prove independent patients; keep gene universe/normalization fixed; confidence vs prediction interval distinction |
| [SyntheVAEiser 2024](https://link.springer.com/article/10.1186/s13059-024-03431-3) | Gene-expression VAE; leave-cancer-out pretraining then fine-tuning on 40 real samples of the target cancer | Few-shot adaptation, not zero-shot transfer to untouched PBTP; not methylation-HRD evidence |
| [GDC methylation pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_Pipeline/) | Current SeSAMe beta processing and raw/masked intensity concepts | Historical and current beta tracks require provenance/bridge |
| [GDC PanImmune](https://gdc.cancer.gov/about-data/publications/panimmune) | Exact filenames/UUID manifests; HRD/components and matrix | Byte-verified data acquired; upstream HRD caller requires further original-source audit |
| [GDC Cell-of-Origin](https://gdc.cancer.gov/about-data/publications/PanCan-CellOfOrigin) | Public SNP6, ABSOLUTE, 450K and RNA v2 manifest | Fast derived inputs; tissue identity is a confounding concern |

Unavailable detailed sources: 2026 Epi-HRD ScienceDirect page, Precious2GPT publisher page, Chen supplementary regression methods, and missing local references. No unverified performance/hyperparameters are assigned to them. MethylNet is an architecture candidate from prior plans, not an independently reproduced KIDS26 baseline. Verify its primary methods before any implementation.
