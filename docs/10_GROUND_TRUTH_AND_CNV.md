# Ground truth, methylation CNV and allele experiments

## Canonical quantities and information limits

**Documented:** Abkevich (local L01 pp. 3,6) counts LOH segments >15 Mb but shorter than a whole chromosome. LOH includes copy-neutral LOH; total CN=2 cannot distinguish A/B=1/1 from 2/0. Whole-chromosome LOH is excluded. Zero-copy deletion is not evidence of an observable retained allele; delegate edge rules to a pinned canonical implementation.

[Popova's original LST definition](https://aacrjournals.org/cancerres/article-abstract/72/21/5454/576090) counts breaks between adjacent >=10 Mb regions after filtering/smoothing <3 Mb variation, separately by arm and excluding centromeric breaks. [Birkbak's TAI definition](https://pmc.ncbi.nlm.nih.gov/articles/PMC3806629/) concerns allelic imbalance reaching a telomere without crossing the centromere. Telomeric total-CN events are not equivalent.

[scarHRD source](https://github.com/sztup/scarHRD) accepts sequenza-style allelic data or sample/chromosome/start/end/total/A/B segmentation; `scar_score` exposes genome build and LOH size cutoff and invokes preprocessing plus separate HRD/AI/LST functions. Its [entry point](https://raw.githubusercontent.com/sztup/scarHRD/master/R/scar_score.R) returns the unadjusted sum in the reviewed revision. Do not silently label a ploidy-adjusted score the same endpoint. Pin a commit, source files, reference chromosome tables and parameters before canonical recalculation; do not rewrite it from a prose definition alone.

**Inference:** canonical LST may include boundaries between allelic states with unchanged total CN, and array resolution/segmentation further changes event counts. Our methyl_LST_like is a total-CN structural proxy until paired concordance is demonstrated. Counts of sufficiently long adjacent total-CN segments are useful even when they are not canonical LST.

## Target hierarchy

| Rank | TCGA | PBTP | Use |
|---|---|---|---|
| Preferred operational target | Versioned published PanImmune `TCGA.HRD_withSampleID.txt`; actual headers sampleID,ai1,lst1,hrd-loh,HRD; checksum retained | Independent paired tumor/normal WGS allele-specific CN -> pinned scarHRD on matching specimen, with reviewed purity/ploidy | Continuous unadjusted sum and components |
| Secondary reference | Rempel ASCAT-derived SNP6 score resource; ABSOLUTE/ASCAT recalculation subset | WES allele-specific CN/scarHRD with sufficient breadth/QC | Paired agreement; stratify by assay; no automatic pooling |
| Orthogonal HR biology | BRCA1/2/PALB2/RAD51C biallelic loss, WGS HRDetect/CHORD, promoter silencing with dosage | Same; RAD51 functional assay if feasible | Corroboration, not numeric replacement |
| Exploratory | RNA scores, total-CN instability, generic HRR mutations | Model-lineage changes, response associations | Hypothesis generation only |

PanImmune is a convenient existing target, not an error-free biological truth. Confirm original upstream caller/coordinate documentation before describing it specifically as ABSOLUTE- or ASCAT-derived; accompanying ABSOLUTE files alone do not prove its provenance. Rempel explicitly uses ASCAT pan-cancer and Sequenza for the OV WES subset (local L10 pp.10-11). Scores shared through GDC, PanCanAtlas and PanImmune are not automatically independent resources or replicate measurements. Compare overlapping sample-level components, median bias, MAE, concordance, Bland-Altman limits, purity/ploidy effects and missingness. Preserve all source IDs; resolve incompatible labels rather than averaging them.

The [scarHRD paper's paired SNP6/WES comparison](https://discovery.ucl.ac.uk/10053224/1/Migrating%20the%20SNP%20array-based%20homologous%20recombination%20deficiency%20measures%20to%20next%20generation%20sequencing%20data%20of%20breast%20cancer.pdf) reports component correlations 0.73-0.84 and sum 0.87. This supports feasibility, not interchangeability or WES-WGS equivalence in pHGG. WES has uneven target coverage, fewer informative heterozygous loci and less precise breakpoints; WGS generally provides denser evidence, but purity and caller ambiguity remain. Use paired WES/WGS for a bridge; with no bridge, report assay-stratified results and acknowledge that prediction error includes reference error. Do not recalibrate on locked PBTP.

## CNV caller and benchmark

Choose conumee2 with fixed assay-compatible normal controls and pinned manifest/build. [Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10868300/) reports LUSC n=367 SNP-array comparison, correlation 0.91; focal deletion performance is weaker than broad agreement (LGG n=239). It uses tangent normalization and weighted segmentation. This is evidence for total-CN feasibility, not HRDsum recovery. Comparison caller: SeSAMe CNV on identical raw samples, only if core pipeline is stable. CopyNumber450kCancer is a baseline-correction alternative with manual-review caveats (local L42 pp.1-2).

Inputs: IDAT red/green pairs, compatible normals, exact array generation, manifest/blacklist, genome build and sample provenance. Beta=M/(M+U) loses absolute M+U. Do not back-calculate total intensity. Low-level M/U is sufficient only if the caller's preprocessing/reference assumptions are met.

Benchmark 20-30 paired public samples if available, selected using assay metadata across purity/ploidy and broad-CN burden; reserve a validation subset before tuning caller settings. Tune reference construction without outcome HRD. Compare on intersected callable genomic intervals, length-weighted agreement, robust log-ratio slope/bias, arm-event sensitivity/specificity, segment/breakpoint precision-recall at predefined 0.5/1/5 Mb tolerances, FGA, event-size strata and LST-like vs genomic LST count differences. Do not count genomic bins as independent tumors. Patient bootstrap; per-patient plots. Homologous chromosomes/centromeres/build exclusions must be consistent.

Initial features: covered-base FGA, gain/loss fraction, amplitude-weighted burden, long-event count, segment-size distribution and arm-specific transition count. Store denominator/callable fraction. Add focal amplification/telomeric proxies only with sufficient resolution. Genome-wide SNP6 correlation alone does not certify breakpoints, allele calls or focal events. Predeclare pragmatic pilot gates (e.g. median broad-CN correlation >=0.8 with acceptable visual QC); choose thresholds from development, never PBTP outcome. If unmet, drop CN predictor and report benchmark failure.

## EPIC allelic subproject: explicit experimental status

[Zhou, Laird and Shen's primary probe-characterization paper](https://zwdzwd.s3.amazonaws.com/papers/2017NAR.pdf) describes explicit SNP probes and additional Infinium-I color-switch polymorphisms. The primary paper (pp.1-2) confirms 59 explicit plus 1052 additional loci for original EPIC. Transfer of those counts and probe identities to the deployed EPIC version requires its exact manifest and raw channel availability. [Author annotation resource](https://zwdzwd.github.io/InfiniumAnnotation) provides SNP/strand/channel annotations; do not use genotype-sensitive probes as ordinary methylation predictors.

**Recommendation:** use matched normal WGS/SNP genotypes to identify genuinely heterozygous loci. A homozygous tumor alone cannot distinguish germline homozygosity from acquired LOH. Retain raw in-band/out-of-band channels before masking, allele-orient per manifest, calibrate A/B ratios on normals, and assess probe reliability versus tumor purity, CN, genotype and tissue methylation. Explicit SNP and color-switch panels should be evaluated separately; avoid misreading bisulfite conversion or allele-specific methylation as allelic CN.

Blindly benchmark locus dosage/BAF concordance, informative loci per chromosome arm, long LOH (>15 Mb) sensitivity/specificity, copy-neutral LOH, telomeric AI and false positives in balanced gains. Compare against WGS allelic segments. Split patients before probe selection/calibration; no PBTP locked-test development. Require adequate informative-locus spacing and held-out event sensitivity/specificity, predeclared for the intended research use, before further investment. Sparse genotyping/ancestry success is not proof of genome-wide allele-specific tumor reconstruction. Failure is a valuable stop result.

## Contradictions preserved

- Ball (L14 pp.3-4) names OTR total-CN-derived measures TAI/LST and bridges them to HRDest. This does not remove TAI's allelic definition. Treat as empirical surrogate evidence, not mechanistic equivalence.
- Chen (L37 pp.5-6) uses MLPA probability >0.5 and differing cohort preprocessing; Rempel (L10) uses genomic sum and context-specific thresholds. Do not collapse their endpoints.
- Commercial GII>0 (L41), percent LOH, GIS>=42 and CHORD/HRDetect probabilities are different scales.
- HRD scars and restored repair after promoter demethylation/reversion can disagree without either assay being technically wrong (L33/L34/L09).
- A persistent Europe PMC 500 does not by itself establish an embargo; the previous project memory's causal interpretation remains unverified. We did not automate institutional access.

## Differential primary-source update — 2026-09-10

- Birkbak L02 derives TAI from allele-specific SNP-array states, defines it as allelic imbalance reaching a subtelomere without crossing the centromere, allows copy-neutral events, uses probe-density thresholds by platform, and applies ASCAT purity/ploidy correction with a reported 36% purity exclusion. This directly reinforces the information limit of total CN.
- Telli L04 reports operational component rules: subchromosomal LOH greater than 15 Mb; TAI to a subtelomere without centromere crossing; LST transitions between segments greater than 10 Mb after filtering segments under 3 Mb. These rules define a genomic reference implementation, not an automatic conumee output.
- Yuan L12 trains a linear SVM GSS from 28 length/location/type/breakpoint features derived by Sequenza using depth and B-allele frequency with purity/ploidy. Its high extreme-case AUC does not validate total-CN methylation reconstruction or pan-cancer transfer.
- The AMP/ACCC/CAP consensus L15 requires laboratories to document component definitions/algorithms and tumor-specific thresholds, include BRCA-wild-type HRD-positive validation, and account for reversion and biallelic status. BRCA1/RAD51C methylation may be reported as a mechanism, but biological thresholds remain context-dependent.

Final terminology is therefore fixed: **reference HRDsum** for independently measured allele-aware LOH+LST+TAI; **predicted reference HRDsum** for the CpG model output; **methylation-derived total-CN structural instability score**, **LST-like breakpoint burden**, or the exact feature names for CN surrogates; and **telomeric CN-alteration burden** rather than TAI when allelic imbalance is unavailable.
