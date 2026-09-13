# Differential literature reconciliation — 2026-09-10

| Paper | Previously reviewed? | Newly available? | Methods significance | Potential to change current plan? |
|---|---|---|---|---|
| L02 Birkbak, *Telomeric allelic imbalance…* | No | Yes | Defines TAI from allele-specific SNP arrays/ASCAT and purity gates | High: confirms allelic requirement |
| L03 Telli, *HRD Score Predicts Response…* | No | Yes | Operational HRDsum and threshold-42 validation | High: tightens target provenance |
| L04 Telli, *Intratumor Heterogeneity…* | No | Yes | Exact LOH/LST/TAI rules; multi-biopsy ICC | Medium: adds lineage/specimen sensitivity |
| L07 Macintyre, *Copy number signatures…* | No | Yes | Absolute-CN feature distributions and NMF signatures | Medium: refines secondary CN branch |
| L12 Yuan, *Genomic Scar Score…* | No | Yes | Sequenza allele-specific CN, 28 features, linear SVM | High: confirms genomic reference requirement |
| L15 Moncur et al., HRD laboratory consensus | No | Yes | Target, algorithm, tumor-specific cutoff, reversion, and methylation reporting guidance | High: refines claims and provenance |
| L19 Doig et al., HR repair deficiency overview | No | Yes | Clinical/pathology synthesis of non-equivalent HRD assays | Medium: confirms cautious naming |
| L23 PARP-inhibitor biomarker meta-analysis | No | Yes | Outcome heterogeneity across HRD biomarker definitions | Medium: confirms no universal clinical cutoff |
| L24 Lord & Ashworth, *BRCAness revisited* | No | Yes | Biological framework; scars versus current function | Medium: confirms scope |
| L27 PARPi/chemotherapy biomarker meta-analysis | No | Yes | Treatment-level evidence synthesis | Low: no methylation implementation |
| L28 ARIEL3 | No | Yes | BRCA/LOH clinical strata and assay-specific endpoint | Medium: target context only |
| L31 Ruscito et al., BRCA1 promoter methylation | No | Yes | Locus-specific methylation in HGSOC | Low: annotation, not HRDsum model |
| L36 HOXA9 methylated ctDNA | No | Yes | Serial locus-specific ddPCR in 32 BRCA-mutated ovarian cancers | Low: monitoring endpoint |
| L39 BRCA1 promoter methylation testing | No | Yes | Direct comparison with myChoice GIS | Medium: promoter annotation, not replacement target |
| L44 Capper et al., CNS methylation classifier | No | Yes | 2,801-case reference, explicit 450K/EPIC mask, calibration, prospective validation | High: refines bridge and n-of-1 design |
| L46 Koelsche et al., methylation classifiers | No | Yes | Shared-platform/QC filtering precedent | Medium: confirms bridge design |
| L47 pediatric brain-tumor methylation review | No | Yes | Pediatric diagnostic context | Low: no HRD method |
| L48 Pajtler et al., ependymal classification | No | Yes | 500 tumors; stable lineage structure | Medium: reinforces cancer-aware validation |
| L49 Sturm et al., H3F3A/IDH1 glioblastoma | No | Yes | Explicit probe mask; strong glioma lineage signal | Medium: reinforces transfer controls |
| L52 Kuschel et al., nanopore classifier | No | Yes | Shared-CpG coverage gate and unclassifiable outcomes | High: refines abstention requirements |
| L59 2026 synthetic-cancer-data review | No | Yes | Separates fidelity, utility, and privacy validation | Medium: confirms deferral |
| L61 omicsGAN AML | No | Yes | 173-case multi-omics WGAN; internal augmentation comparison | Low: no HRD/external transfer |
| L63 synthetic-data oncology review | No | Yes | General opportunities and validation risks | Low: no direct benchmark |
| L65 Krause et al., histology GAN | No | Yes | Real versus synthetic versus combined training on real test images | Low: useful experiment design, wrong modality |
| L70 chordoma HRD/PARP study | No; absent from matrix | Yes | Matched tumor-normal sequencing and rare-tumor HRD scars | Medium: biological motivation only |

## Inventory reconciliation

The prior frozen inventory listed 35 local PDFs. The current folder contains 60 PDFs: the same 35 plus 25 newly available files. Of the 25, 24 resolve previously unavailable rows in the supplied 69-entry map and L70 is additional local evidence absent from that map. SHA-256 comparison found no exact duplicate PDFs. No newly added file is a supplementary-only appendix; several publisher/NIH manuscript filenames are opaque or truncated main-article names, so the matrix records citations separately from local paths.

The machine-readable extraction for every new paper is in [literature_evidence.csv](literature_evidence.csv): cohort, N, assay, input, target, component definitions, allele-specific method, purity/ploidy, preprocessing, model, feature selection, training, validation, transfer, metrics, uncertainty, single-sample use, code, synthetic method, limitations, and KIDS26 consequence. Unreported details are marked absent or not established rather than inferred.

## What Changed Because of the Newly Added Literature?

### 1. HRD target provenance became a release-blocking requirement

**Previous conclusion:** predict an independently measured continuous HRDsum and avoid treating all commercial or literature HRD measures as interchangeable.

**New evidence:** Telli's score is an unweighted LOH+TAI+LST sum, with the threshold of 42 derived in an independent BRCA-deficient reference set; the new laboratory consensus asks reports to state each component, algorithm, locus distribution, specimen context, and tumor-specific cutoff.

**Updated conclusion:** the architecture is unchanged, but a label-source/caller/build/parameter record is mandatory. A pan-cancer threshold of 42 is exploratory unless validated for that exact assay and cancer context.

**Implementation consequence:** model the continuous published reference first; report components and provenance; do not expose a portable clinical HRD-high probability without independent calibration.

### 2. The total-CN limitation is now more strongly supported

**Previous conclusion:** ordinary methylation-derived total CN cannot recover canonical LOH or TAI; LST-like information is more plausible.

**New evidence:** Birkbak defines telomeric allelic imbalance using unequal allele copy number, allows copy-neutral events, and uses ASCAT purity/ploidy correction. Telli's intratumor paper defines LOH as subchromosomal events over 15 Mb, TAI as allelic imbalance reaching a subtelomere without crossing the centromere, and LST as breaks between segments over 10 Mb after filtering segments under 3 Mb. Yuan's GSS also starts from Sequenza depth plus B-allele frequency.

**Updated conclusion:** unchanged. Theoretical genotype signal from explicit SNP probes or Type-I color switching is not a validated genome-wide canonical-scar implementation in the supplied literature.

**Implementation consequence:** call any methylation-total-CN output a **methylation-derived total-CN structural instability score**, **LST-like breakpoint burden**, or named feature set. Reserve *HRDsum* for the allele-aware reference or a clearly labeled prediction of it.

### 3. The feature bridge and n-of-1 contract are stricter

**Previous conclusion:** freeze a shared 450K/EPIC feature list and require OOD handling.

**New evidence:** Capper's diagnostic system uses per-sample preprocessing, an explicit cross-platform/QC mask, calibrated scores, and prospective validation. Kuschel requires shared-CpG coverage and permits unclassifiable results. Glioma/ependymoma studies demonstrate large lineage signals.

**Updated conclusion:** use an exact versioned manifest intersection, record hashes and coverage, and predeclare abstention. Cancer-aware evaluation remains indispensable.

**Implementation consequence:** the repository now contains a deterministic bridge builder and an EPIC-v1 candidate mask; the PBTP EPIC generation must be confirmed before freeze.

### 4. The CN branch is better specified, not promoted

**Previous conclusion:** gate total-CN structural features and fusion behind technical agreement.

**New evidence:** Macintyre derives interpretable absolute-CN distributions and NMF signatures with extensive stability checks, but depends on purity/ploidy-aware genomic CN and warns that resolution/preprocessing matter.

**Updated conclusion:** the branch remains secondary. Start with simple prespecified burden/breakpoint features; avoid a new NMF discovery exercise during the event.

**Implementation consequence:** compare CpG-only, total-CN-feature-only, and fusion on identical paired patients and splits; promote only by predeclared development improvement.

### 5. Synthetic augmentation remains deferred

**Previous conclusion:** establish the real-data baseline and transfer problem before augmentation.

**New evidence:** the 2026 review requires separate fidelity, utility, and privacy evidence. omicsGAN reports an internal AML gain in 173 cases but lacks an HRD target, external cohort, or 450K-to-EPIC test. The histology GAN supports testing augmentation only on untouched real data, but its modality does not transfer.

**Updated conclusion:** unchanged. Synthetic data would mostly add development risk during the three-day event.

**Implementation consequence:** at most, run one predeclared train-only PCA-space interpolation sensitivity after the minimum and strong results exist. Compare real-only versus real-plus-synthetic training on identical untouched real validation; never claim synthetic samples as independent validation.

## Central HRDsum assessment

### HRD-LOH

**Documented fact:** canonical HRD-LOH requires allele-specific states. The new papers use SNP/sequence B-allele information and, where specified, ASCAT or Sequenza with purity/ploidy correction. Copy-neutral LOH is invisible in balanced diploid total CN.

**Inference:** 450K/EPIC contains some explicit SNP and genotype-sensitive intensity information, so an experimental BAF-like assay is theoretically possible with raw channels, manifest-aware allele orientation, informative heterozygous loci, and preferably matched normal genotypes. No supplied paper validates this as a canonical genome-wide HRD-LOH implementation.

**Recommendation:** do not calculate or name canonical HRD-LOH from conumee/total-CN segments. Keep the raw-channel allelic project separate and benchmark it blindly against WGS/SNP-array allele-specific segments before any HRD claim.

### LST

**Documented fact:** canonical LST counts transitions between adjacent large segments after small-segment filtering; the newly reviewed implementation uses segments over 10 Mb after removing segments under 3 Mb. Total-CN methylation segmentation can expose large breakpoints, but caller smoothing, centromeres, purity, ploidy, probe density, and breakpoint tolerance alter counts.

**Inference:** a methylation **LST-like breakpoint burden** is plausible if broad-CN concordance and breakpoint recovery are demonstrated on paired specimens. It is not automatically numerically equivalent to a canonical allele-aware implementation.

**Recommendation:** freeze segment smoothing, minimum 3 Mb filtering, 10 Mb large-segment rule, centromere/build handling, and breakpoint tolerances before validation; report agreement and bias against the genomic reference.

### TAI

**Documented fact:** TAI requires allelic imbalance extending to the subtelomere without crossing the centromere; copy-neutral and copy-number-altered events can qualify. A telomeric total-CN alteration alone does not establish allelic imbalance.

**Recommendation:** do not call telomeric total-CN burden TAI. Use **telomeric CN-alteration burden** unless an allele-aware raw-channel method passes external event-level validation.

## Model-family decision

| Candidate | Decision | Reason |
|---|---|---|
| Elastic net on shared CpGs | Primary | Fits high-p/medium-n data, sparse/regularized, fast nested tuning, interpretable, and easy to freeze for one sample |
| Elastic net on prespecified total-CN features | Secondary #1 | Mechanistically interpretable structural comparison after technical gate |
| Elastic-net fusion | Secondary #2 | Clean incremental test on the same paired samples/splits |
| Ridge/lasso | Included in the elastic-net alpha grid | Useful endpoints, no separate pipeline needed |
| RF/boosting/SVM/PLS/PCA regression | Deferred | More tuning or preprocessing burden; no new HRD transfer evidence justifies displacement |
| Neural nets/autoencoders/VAE/GAN/diffusion | Post-hackathon | Sample size, leakage, cross-platform, and validation burdens exceed the evidence |
| Pathway/ssGSEA-like scores | Exploratory interpretation | May summarize biology but no supplied methylation-HRDsum benchmark supports primary use |
| Stacking/ensembles | Post-hackathon | Requires extra out-of-fold infrastructure and a stable base-model benchmark |

**Primary hackathon model:** shared-probe CpG elastic-net regression to continuous independently measured reference HRDsum.

**Secondary model #1:** methylation-derived total-CN structural features to reference HRDsum, after paired technical validation.

**Secondary model #2:** CpG plus total-CN feature fusion on identical paired cases and splits.

**Deferred post-hackathon methods:** component models, nonlinear learners, pathway models, ensembles, raw-channel experimental BAF, and synthetic/deep generative augmentation.

The newly added literature provides no sufficient reason to replace elastic net.
