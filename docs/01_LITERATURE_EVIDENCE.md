# Literature evidence: paper-by-paper digest

Evidence labels distinguish source claims from KIDS26 judgments. Unknown is not a negative finding. Page numbers are physical PDF pages including covers. Original CSV titles/years are retained, with discrepancies noted. This is a methods-focused first pass, not a systematic review or a complete supplementary-methods audit.

## L01 - Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer.

Source: ../papers/2012_Patterns_of_genomic_loss_of_heterozygosity_predict_homologous_recombination_repa.pdf; pp. 1-4,6

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** The original HRD score counts LOH regions; it is not the later three-component sum. Scars can persist after reversion.

**Inference:** HRD score is an ambiguous historical name.

**Recommendation:** Store component names and assay provenance explicitly.

## L02 - Telomeric allelic imbalance indicates defective DNA repair and sensitivity to DNA-damaging agents.

Source: ../papers/nihms-513878.pdf; pp. 1-8; Methods pp. 6-8

Status: Newly available author manuscript; consequential methods reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** TAI was derived from SNP-array allele signals with ASCAT, not total CN. Copy-neutral LOH is explicitly a possible source of telomeric AI.

**Inference:** Neither ordinary beta values nor methylation-derived total CN contains the demonstrated input needed for canonical TAI.

**Recommendation:** Retain canonical TAI only for allele-aware reference data; call total-CN telomeric features telomeric CN alteration burden.

## L03 - Homologous Recombination Deficiency (HRD) Score Predicts Response to Platinum-Containing Neoadjuvant Chemotherapy in Patients with Triple-Negative Breast Cancer.

Source: ../papers/nihms-1052067.pdf; pp. 2-11; Methods pp. 4-7

Status: Newly available author manuscript; main methods reviewed, component supplement unavailable

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** The canonical combined score is an unweighted allele-aware LOH+TAI+LST sum. The clinical cutoff was trained in breast/ovarian BRCA-deficient tumors, not pediatric CNS tumors.

**Inference:** A value of 42 is not automatically portable across cancers or assays.

**Recommendation:** Predict continuous reference HRDsum and do not emit an HRD-high probability without tumor-specific calibration.

## L04 - Intratumor Heterogeneity of Homologous Recombination Deficiency in Primary Breast Cancer.

Source: ../papers/1193.pdf; pp. 1-6

Status: Newly available full article; methods/results reviewed

Delta category: REFINES CURRENT PLAN

**Documented fact:** Allele-aware HRDsum was spatially reproducible in this small breast cohort, but only 31 patients had scores.

**Inference:** Lineage grouping remains mandatory even when scar scores appear stable within one adult tumor type.

**Recommendation:** Keep patient/lineage grouped validation; do not treat related PBTP models as independent.

## L05 - HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures.

Source: ../papers/2017_HRDetect_is_a_predictor_of_BRCA1_and_BRCA2_deficiency_based_on_mutational_signat.pdf; pp. 2,4-7

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** HRDetect integrates several mutation classes and estimates BRCA deficiency probability, not HRDsum.

**Inference:** WES does not recover all WGS rearrangement information.

**Recommendation:** Orthogonal WGS comparator only.

## L06 - A mutational signature reveals alterations underlying defective homologous recombination repair.

Source: https://pubmed.ncbi.nlm.nih.gov/?term=A+mutational+signature+reveals+alterations+underlying+defective+homologous+recombination+repair;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L07 - Copy number signatures and mutational processes in cancer.

Source: ../papers/emss-78219.pdf; pp. 1-12; figures/supplement references thereafter

Status: Newly available author manuscript; main methods/results reviewed

Delta category: REFINES CURRENT PLAN

**Documented fact:** The seven signatures were derived from total/absolute CN features with strong quality filtering; the paper explicitly warns that preprocessing and ploidy alter signatures.

**Inference:** This supports a total-CN structural branch but cannot recover LOH or TAI without alleles.

**Recommendation:** Keep total-CN features secondary and name them as structural/CN signatures, not HRDsum components.

## L08 - Signatures of copy number alterations in human cancer.

Source: ../papers/2022_Signatures_of_copy_number_alterations_in_human_cancer.pdf; pp. 1-2,9

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Describes 21 CN signatures; 17 attributed to processes including doubling, LOH and HRD. Feature construction includes allelic information.

**Inference:** Total-intensity CN cannot automatically reproduce allele-specific signature channels.

**Recommendation:** Borrow burden/length concepts; label partial signatures as surrogates.

## L09 - Pan-cancer landscape of homologous recombination deficiency.

Source: ../papers/2020_Pan_cancer_landscape_of_homologous_recombination_deficiency.pdf; pp. 1-3,9

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** CHORD validates across cancers and discusses historical scars, mutation-count QC and variant-caller dependence.

**Inference:** QC-failed low-mutation pediatric genomes cannot be called HR-proficient.

**Recommendation:** Keep probability and QC as separate orthogonal evidence.

## L10 - Pan-cancer analysis of genomic scar patterns caused by homologous repair deficiency (HRD).

Source: ../papers/2022_Pan_cancer_analysis_of_genomic_scar_patterns_caused_by_homologous_repair_deficie.pdf; pp. 1,6-8,10-11

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Pan-cancer scores use SNP-array ASCAT; WES is an ovarian comparison. Cutpoints and purity effects differ by cancer.

**Inference:** Prior reading-map wording about WES pan-cancer labels is misleading.

**Recommendation:** One versioned target; paired comparison before pooling labels.

## L11 - Inferring Homologous Recombination Deficiency of Ovarian Cancer From the Landscape of Copy Number Variation at Subchromosomal and Genetic Resolutions.

Source: ../papers/2022_Inferring_Homologous_Recombination_Deficiency_of_Ovarian_Cancer_From_the_Landsca.pdf; pp. 1-3

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** 8q24.2/5q13.2/19q12 events associate with OV HRD. External labels differ. PDF date is December 2021; CSV says 2022.

**Inference:** Cross-cohort association is not calibration of the same continuous endpoint.

**Recommendation:** Keep ovarian loci exploratory.

## L12 - Genomic Scar Score: A robust model predicting homologous recombination deficiency based on genomic instability.

Source: ../papers/BJOG - 2022 - Yuan - Genomic Scar Score  A robust model predicting homologous recombination deficiency based on genomic.pdf; pp. 1-7; Methods pp. 2-4

Status: Newly available full article; methods/results reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** Despite broad wording about chromosomal CN, the model requires Sequenza BAF-derived allele-specific CN and purity/ploidy.

**Inference:** Its high AUC does not show canonical scars can be obtained from methylation total CN.

**Recommendation:** Do not replace elastic-net CpG regression; consider only as an allele-aware genomic comparator.

## L13 - Clinical evaluation of a low-coverage whole-genome test for detecting homologous recombination deficiency in ovarian cancer.

Source: ../papers/2024_Clinical_evaluation_of_a_low_coverage_whole_genome_test_for_detecting_homologous.pdf; pp. 1,3; local preprint

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** The local file is a medRxiv preprint version of a shallow-WGS HRD evaluation.

**Inference:** Its clinical agreement cannot establish methylation-CN equivalence.

**Recommendation:** Use paired-assay validation design; verify final version before citing clinical conclusions.

## L14 - Leveraging Off-Target Reads in Panel Sequencing for Homologous Recombination Repair Deficiency Screening in Tumor.

Source: ../papers/2024_Leveraging_Off_Target_Reads_in_Panel_Sequencing_for_Homologous_Recombination_Rep.pdf; pp. 2-4,6

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Calls OTR features TAI+LST despite sparse SNP information; multiplies by 1.2 for HRDest (1.3 in TCGA/WES-only comparison). Reference cases span a selected wide score range.

**Inference:** TAI terminology does not prove allele recovery. Range selection and calibration/evaluation reuse limit transportability.

**Recommendation:** Name our total-CN telomere feature telomeric_CNV_burden; require paired allele-specific validation.

## L15 - Recommendations for Clinical Molecular Laboratories for Detection of Homologous Recombination Deficiency in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, Association of Cancer Care Centers, and College of American Pathologists.

Source: ../papers/1-s2.0-S1525157825001369-main.pdf; pp. 1-16; recommendations pp. 10-14

Status: Newly available consensus article; recommendations and validation sections reviewed

Delta category: REFINES CURRENT PLAN

**Documented fact:** The consensus recommends declaring specimen/tumor-content requirements, validating BRCA-wild-type HRD-positive samples, documenting each scar and its detection method, and using tumor-specific cutoffs. BRCA1/RAD51C promoter methylation may be considered, but relevant thresholds remain unclear.

**Inference:** A pan-cancer binary threshold or unlabeled composite would be scientifically weak for KIDS26.

**Recommendation:** Keep continuous output primary; add assay provenance and tumor-specific cutoff warnings; make promoter methylation an orthogonal annotation.

## L16 - Patient Assessment and Therapy Planning Based on Homologous Recombination Repair Deficiency.

Source: ../papers/2023_Patient_Assessment_and_Therapy_Planning_Based_on_Homologous_Recombination_Repair.pdf; pp. 1-3

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Reviews non-equivalent HRD assays and breast/ovarian origin of GIS 42.

**Inference:** Public scar sums cannot inherit commercial clinical claims.

**Recommendation:** Separate phenotype, assay, cutoff and clinical endpoint.

## L17 - Biomarkers for Homologous Recombination Deficiency in Cancer.

Source: https://pubmed.ncbi.nlm.nih.gov/29788099/;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L18 - Biomarkers for Homologous Recombination Deficiency in Cancer.

Source: ../papers/2021_Biomarkers_for_Homologous_Recombination_Deficiency_in_Cancer.pdf; pp. 1,8

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Discusses combined scars and limitations for current repair status.

**Inference:** HRDsum is a historical phenotype.

**Recommendation:** Add functional corroboration later.

## L19 - Homologous Recombination Repair Deficiency: An Overview for Pathologists.

Source: ../papers/1-s2.0-S089339522200480X-main.pdf; pp. 1-7

Status: Newly available review; relevant definitions/limitations reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** Low tumor purity and heterogeneity impair mutation, scar and signature detection; scars can reflect etiology without proving current HR function.

**Inference:** KIDS26 requires purity-aware evaluation and cannot make therapy claims from predicted historical scars.

**Recommendation:** Preserve research-only continuous phenotype language and orthogonal functional corroboration roadmap.

## L20 - Homologous recombination deficiency in breast cancer: Implications for risk, cancer development, and therapy.

Source: ../papers/2021_Homologous_recombination_deficiency_in_breast_cancer_Implications_for_risk_cance.pdf; pp. 13,18; manuscript

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Discusses breast HRD and persistent scars after repair-restoring reversion.

**Inference:** Therapy response and scar prediction differ.

**Recommendation:** No therapeutic eligibility from prototype.

## L21 - ESMO recommendations on predictive biomarker testing for homologous recombination deficiency and PARP inhibitor benefit in ovarian cancer.

Source: https://pubmed.ncbi.nlm.nih.gov/33004253/;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L22 - Ovarian Cancer Therapy: Homologous Recombination Deficiency as a Predictive Biomarker of Response to PARP Inhibitors.

Source: ../papers/2022_Ovarian_Cancer_Therapy_Homologous_Recombination_Deficiency_as_a_Predictive_Bioma.pdf; pp. 1-6

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Summarizes ovarian PARPi trials and HRD assays at publication.

**Inference:** Historical drug indications are time dependent.

**Recommendation:** Biological background only.

## L23 - Comparison of PARPi efficacy according to homologous recombination deficiency biomarkers in patients with ovarian cancer: a systematic review and meta-analysis.

Source: ../papers/113706-PB12-2453-R2.pdf; pp. 1-7

Status: Newly available systematic review/meta-analysis; methods/results reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** The meta-analysis treats Foundation genomic LOH and myChoice LOH+TAI+LST as different biomarker definitions.

**Inference:** Correlated scar assays are not interchangeable labels.

**Recommendation:** Keep one versioned reference endpoint and do not inherit clinical claims.

## L24 - BRCAness revisited.

Source: ../papers/nrc.2015.21.pdf; pp. 1-10

Status: Newly available review; relevant biology/scar sections reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** Genomic scars portray tumor history and may be collateral consequences rather than current drivers.

**Inference:** A scar predictor is not a direct functional-HR measurement.

**Recommendation:** Maintain historical-scar language and separate current resistance mechanisms.

## L25 - PARP inhibitors: Synthetic lethality in the clinic.

Source: ../papers/2017_PARP_inhibitors_Synthetic_lethality_in_the_clinic.pdf; pp. 1-4

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Reviews PARP synthetic lethality, unrelated to synthetic-data generation.

**Inference:** Mechanism motivates but does not validate prediction.

**Recommendation:** Background only.

## L26 - PARP Inhibitors as a Therapeutic Agent for Homologous Recombination Deficiency in Breast Cancers.

Source: ../papers/2019_PARP_Inhibitors_as_a_Therapeutic_Agent_for_Homologous_Recombination_Deficiency_in_.pdf; pp. 1,10-11

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Reviews breast PARPi mechanisms and trials.

**Inference:** Adult efficacy does not validate pediatric scar predictions.

**Recommendation:** Research phenotype claims only.

## L27 - DNA repair biomarkers to guide usage of combined PARP inhibitors and chemotherapy: A meta-analysis and systematic review.

Source: ../papers/1-s2.0-S1043661823002839-main.pdf; pp. 1-11

Status: Newly available systematic review/meta-analysis; methods/results reviewed

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** The clinical meta-analysis pools distinct DNA-repair biomarkers and does not validate a common numeric HRD scale.

**Inference:** Clinical enrichment cannot validate predicted TCGA HRDsum.

**Recommendation:** Do not use treatment response as the hackathon training target.

## L28 - Rucaparib maintenance treatment for recurrent ovarian carcinoma after response to platinum therapy (ARIEL3): a randomised, double-blind, placebo-controlled, phase 3 trial.

Source: ../papers/nihms956813.pdf; pp. 1-15

Status: Newly available trial manuscript; biomarker design/results reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** ARIEL3 used genome-wide LOH, not the three-component HRDsum.

**Inference:** Its clinical threshold cannot calibrate KIDS26.

**Recommendation:** Preserve endpoint names and provenance.

## L29 - Homologous Recombination Repair Gene Mutations to Predict Olaparib Plus Bevacizumab Efficacy in the First-Line Ovarian Cancer PAOLA-1/ENGOT-ov25 Trial.

Source: ../papers/2023_Homologous_Recombination_Repair_Gene_Mutations_to_Predict_Olaparib_Plus_Bevacizu.pdf; pp. 1,3,8

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Non-BRCA HRR mutation panels did not predict maintenance benefit; biallelic frequencies varied.

**Inference:** Any HRR mutation is not a universal positive label.

**Recommendation:** Use curated biallelic defects as corroboration.

## L30 - Clinical and molecular characteristics of ARIEL3 patients who derived exceptional benefit from rucaparib maintenance treatment for high-grade ovarian carcinoma.

Source: https://pubmed.ncbi.nlm.nih.gov/36273926/;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L31 - BRCA1 gene promoter methylation status in high-grade serous ovarian cancer patients--a study of the tumour Bank ovarian cancer (TOC) and ovarian cancer diagnosis consortium (OVCAD).

Source: ../papers/1-s2.0-S0959804914006492-main.pdf; pp. 1-7

Status: Newly available full article; assay and outcome sections reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** BRCA1 promoter methylation alone was not associated with PFS or OS in this cohort.

**Inference:** Single-gene methylation is neither a universal HRD label nor a substitute for HRDsum.

**Recommendation:** Keep promoter status as an orthogonal feature/interpretation, not primary ground truth.

## L32 - BRCA1 promoter methylation is a marker of better response to platinum-taxane-based therapy in sporadic epithelial ovarian cancer.

Source: ../papers/2014_BRCA1_promoter_methylation_is_a_marker_of_better_response_to_platinum_taxane_bas.pdf; pp. 1-3

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** BRCA1 methylation investigated in two ovarian cohorts receiving platinum/taxane.

**Inference:** Promoter positivity is assay-specific.

**Recommendation:** Mechanistic context only.

## L33 - Methylation of all BRCA1 copies predicts response to the PARP inhibitor rucaparib in ovarian carcinoma.

Source: ../papers/2018_Methylation_of_all_BRCA1_copies_predicts_response_to_the_PARP_inhibitor_rucapari.pdf; pp. 7,11-13

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Methylation of all BRCA1 copies and current biopsy status matter; low purity can prevent zygosity inference.

**Inference:** A promoter beta cutoff ignores dosage/purity.

**Recommendation:** Do not use BRCA1 beta positivity as ground truth.

## L34 - Acquired RAD51C Promoter Methylation Loss Causes PARP Inhibitor Resistance in High-Grade Serous Ovarian Carcinoma.

Source: ../papers/2021_Acquired_RAD51C_Promoter_Methylation_Loss_Causes_PARP_Inhibitor_Resistance_in_Hi.pdf; pp. 3,24-25

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** RAD51C methylation loss accompanied restored expression and resistance; scars reflected HRD history.

**Inference:** Derived models can evolve repair phenotypes.

**Recommendation:** Group lineages; prioritize patient tumors for external testing.

## L35 - High-level tumour methylation of BRCA1 and RAD51C is required for homologous recombination deficiency in solid cancers.

Source: ../papers/2024_High_level_tumour_methylation_of_BRCA1_and_RAD51C_is_required_for_homologous_rec.pdf; pp. 1-3

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** High-level BRCA1/RAD51C methylation with LOH relates to HRD; lower-level RAD51C methylation in LGG/stomach is often associated with CIMP instead.

**Inference:** CNS hypermethylation can confound adult HRD predictors.

**Recommendation:** Stratify errors by subtype and purity; do not equate hypermethylation with HRD.

## L36 - Circulating HOXA9-methylated tumour DNA: A novel biomarker of response to poly (ADP-ribose) polymerase inhibition in BRCA-mutated epithelial ovarian cancer.

Source: ../papers/1-s2.0-S0959804919308317-main.pdf; Methods pp. 3-5; Results pp. 5-8

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** HOXA9 methylated ctDNA was measured by bisulfite ddPCR in 32 BRCA1/2-mutated ovarian cancer patients treated with veliparib.

**Inference:** A treatment-monitoring methylation locus cannot establish that methylation arrays calculate canonical HRDsum.

**Recommendation:** Do not alter the hackathon model family based on this study.

## L37 - Extensive epigenomic dysregulation is a hallmark of homologous recombination deficiency in triple-negative breast cancer.

Source: ../papers/2024_Extensive_epigenomic_dysregulation_is_a_hallmark_of_homologous_recombination_def.pdf; pp. 4-7,10-12

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Predicts MLPA BRCAness probability from methylation, not canonical HRDsum. Two biological cohorts support the association; model details are referred to the supplement.

**Inference:** Supports plausibility, not validated adult-pan-cancer-to-pediatric regression.

**Recommendation:** Obtain supplement before reproducing its tuning; independently implement nested baseline.

## L38 - Development and validation of an epigenetic score for homologous recombination deficiency based on genome-wide DNA methylation profiling in ovarian cancer.

Source: https://www.sciencedirect.com/science/article/pii/S2452014426002013;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L39 - Clinical Relevance of BRCA1 Promoter Methylation Testing in Patients with Ovarian Cancer.

Source: ../papers/3124.pdf; Methods pp. 3-5; Results pp. 5-8

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: REFINES CURRENT PLAN

**Documented fact:** High and low BRCA1 methylation strata were compared with myChoice GIS in ovarian tumors.

**Inference:** Promoter methylation can mark an HRD mechanism but is not a sufficient universal HRDsum surrogate.

**Recommendation:** Report BRCA1/RAD51C promoter features as annotations or secondary features, not the sole target definition.

## L40 - BRCA1 promoter methylation predicts PARPi response in ovarian cancer: insights from the KOMET study.

Source: ../papers/2025_BRCA1_promoter_methylation_predicts_PARPi_response_in_ovarian_cancer_insights_fr.pdf; pp. 1,3,5

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Retrospective study relates BRCA1 methylation to outcome, with missing GIS and treatment selection.

**Inference:** Not validation of a methylation-wide continuous score.

**Recommendation:** Biological context only.

## L41 - A Novel Droplet Digital PCR Assay for BRCA1 and RAD51C Methylation: Advancing Homologous Recombination Deficiency Detection in Ovarian Cancer.

Source: ../papers/2025_A_Novel_Droplet_Digital_PCR_Assay_for_BRCA1_and_RAD51C_Methylation_Advancing_Hom.pdf; p. 1

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Assesses promoter methylation against commercial GII-based HRD.

**Inference:** Commercial scales cannot be pooled numerically.

**Recommendation:** Store assay-specific targets separately.

## L42 - CopyNumber450kCancer: baseline correction for accurate copy number calling from the 450k methylation array.

Source: ../papers/2016_CopyNumber450kCancer_baseline_correction_for_accurate_copy_number_calling_from_the_.pdf; pp. 1-2

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Explicitly notes lack of SNP-array B-allele frequency; this is a baseline correction method.

**Inference:** Label-visible manual correction can leak reference information.

**Recommendation:** Locked baseline/QC; blind to HRD outcomes.

## L43 - MethPed: a DNA methylation classifier tool for the identification of pediatric brain tumor subtypes.

Source: ../papers/2015_MethPed_a_DNA_methylation_classifier_tool_for_the_identification_of_pediatric_br.pdf; pp. 1-2,8

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** MethPed provides saved-forest single-sample classification of pediatric tumor types.

**Inference:** Subtype signal is both useful and a potential HRD confounder; nesting must be independently enforced.

**Recommendation:** Borrow frozen inference, not accuracy expectations.

## L44 - DNA methylation-based classification of central nervous system tumours.

Source: ../papers/nihms-942946.pdf; Methods pp. 18-28; Results pp. 5-12; extended methods/tables pp. 29-45

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: REFINES CURRENT PLAN

**Documented fact:** The classifier removed 11,551 sex-chromosome, 7,998 common-SNP-proximate, 3,965 nonunique hg19, and 32,260 non-EPIC probes, retaining 428,799.

**Inference:** Cross-platform compatibility is feasible with a frozen mask, but classification performance does not prove HRDsum transfer.

**Recommendation:** Adopt deterministic bridge provenance and single-sample preprocessing while retaining cancer-held-out HRD validation.

## L45 - Practical implementation of DNA methylation and copy-number-based CNS tumor diagnostics: the Heidelberg experience.

Source: ../papers/2018_Practical_implementation_of_DNA_methylation_and_copy_number_based_CNS_tumor_diag.pdf; p. 1; diagnostic workflow sections

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Integrates methylation, CN profiles and morphology for CNS diagnostics.

**Inference:** CN corroboration does not validate HRD.

**Recommendation:** Show assay QC and applicability with predictions.

## L46 - Methylation classifiers: Brain tumors, sarcomas, and what's next.

Source: ../papers/Genes Chromosomes   Cancer - 2022 - Koelsche - Methylation classifiers  Brain tumors  sarcomas  and what s next.pdf; Methods overview and cross-platform recommendations throughout

Status: Newly available 2026-09-10; methods-relevant sections reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** The review describes restricting analysis to shared 450K/EPIC probes and removing sex-linked, SNP-proximate, and nonunique probes.

**Inference:** The same platform-compatibility principles are appropriate for KIDS26 with exact manifests.

**Recommendation:** Use a versioned shared autosomal general-mask list and record annotation hashes.

## L47 - Advances in the classification of pediatric brain tumors through DNA methylation profiling: From research tool to frontline diagnostic.

Source: ../papers/Cancer - 2018 - Kumar - Advances in the classification of pediatric brain tumors through DNA methylation profiling  From.pdf; Review sections on molecular targets and pediatric glioma biology

Status: Newly available 2026-09-10; relevance review completed

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** The paper reviews targeted therapies and molecular subgroups in pediatric glioma.

**Inference:** Disease context supplies no evidence for substituting a pediatric classifier for HRD regression.

**Recommendation:** Use only for disease-context framing.

## L48 - Molecular Classification of Ependymal Tumors across All CNS Compartments, Histopathological Grades, and Age Groups.

Source: ../papers/nihms-745680.pdf; Methods pp. 15-20; Results pp. 5-12

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** Consensus clustering of the 10,000 most variable CpGs separated ependymal molecular groups in 500 tumors.

**Inference:** A pooled methylation-HRD model could learn lineage unless cancer-aware validation is enforced.

**Recommendation:** Retain leave-one-cancer-type-out evaluation and lineage-aware PBTP grouping.

## L49 - Hotspot mutations in H3F3A and IDH1 define distinct epigenetic and biological subgroups of glioblastoma.

Source: ../papers/1-s2.0-S1535610812003649-main.pdf; Methods pp. 10-14; Results pp. 3-8

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** The analysis retained 438,370 probes after sex, SNP-proximity, and nonunique filtering, then clustered on 8,000 variable probes.

**Inference:** Lineage methylation can inflate pooled HRD prediction when folds mix cancers.

**Recommendation:** Keep supervised preprocessing inside folds and report cancer-held-out results.

## L50 - Transcriptomic and epigenetic profiling of 'diffuse midline gliomas, H3 K27M-mutant' discriminate two subgroups based on the type of histone H3 mutated and not supratentorial or infratentorial location.

Source: ../papers/2018_Transcriptomic_and_epigenetic_profiling_of_diffuse_midline_gliomas_H3_K27M_mutan.pdf; pp. 1,4-5,11

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** H3.1/H3.3 epigenetic differences can persist in cultures; text/figure differ in duplicate count.

**Inference:** Histone subtype and lineage can dominate methylation distances.

**Recommendation:** Group donor models and assess subtype-specific errors.

## L51 - Ultra-fast deep-learned CNS tumour classification during surgery.

Source: ../papers/2023_Ultra_fast_deep_learned_CNS_tumour_classification_during_surgery.pdf; pp. 1-3

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Sturgeon uses sparse simulation and abstention with real validation.

**Inference:** Millions of masks do not increase independent biological N.

**Recommendation:** Borrow dropout tests/abstention.

## L52 - Robust methylation-based classification of brain tumours using nanopore sequencing.

Source: ../papers/Neuropathology Appl Neurobio - 2022 - Kuschel - Robust methylation‐based classification of brain tumours using nanopore.pdf; Methods pp. 3-7; Results pp. 7-12

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: REFINES CURRENT PLAN

**Documented fact:** The method required at least 1,000 shared CpGs and allowed samples to remain unclassified; specificity was 100% among classifiable validation samples.

**Inference:** A KIDS26 n-of-1 system should expose coverage and abstention rather than force every prediction.

**Recommendation:** Predeclare OOD and feature-coverage diagnostics before PBTP application.

## L53 - Rapid DNA methylation-based classification of pediatric brain tumors from ultrasonic aspirate specimens.

Source: ../papers/2024_Rapid_DNA_methylation_based_classification_of_pediatric_brain_tumors_from_ultrason.pdf; pp. 1,8-9

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Pilot cross-assay study derives platform-specific confidence threshold.

**Inference:** Same-cohort cutoff selection limits external calibration claims.

**Recommendation:** Predeclare PBTP thresholds.

## L54 - Rapid brain tumor classification from sparse epigenomic data.

Source: ../papers/2025_Rapid_brain_tumor_classification_from_sparse_epigenomic_data.pdf; pp. 1-3,7,11,15

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** A simple probabilistic classifier enables rapid sparse methylation classification.

**Inference:** Large feature space does not itself justify neural networks.

**Recommendation:** Regularized baseline first.

## L55 - Prospective, multicenter validation of a platform for rapid molecular profiling of central nervous system tumors.

Source: ../papers/2025_Prospective_multicenter_validation_of_a_platform_for_rapid_molecular_profiling_of_.pdf; pp. 1-2,12

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Large diagnostic validation; Rapid-CNS2 methods include sample-specific forest retraining.

**Inference:** Single-sample has different operational meanings.

**Recommendation:** Require unchanged model weights in KIDS26.

## L56 - Nanopore-based random genomic sampling for intraoperative molecular diagnosis.

Source: ../papers/2025_Nanopore_based_random_genomic_sampling_for_intraoperative_molecular_diagnosis.pdf; pp. 1,11,16

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Combines CN and methylation for intraoperative molecular diagnosis.

**Inference:** Fusion benefit for HRD remains unproven.

**Recommendation:** Evaluate fusion on identical paired samples.

## L57 - TUCAN: Ultra-fast methylation-based classification of pediatric solid tumors and lymphomas.

Source: ../papers/2026_TUCAN_Ultra_fast_methylation_based_classification_of_pediatric_solid_tumors_and.pdf; pp. 2,4,6,13-14,40;preprint v1

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** TUCAN is a preprint with platform intersection and real validation plus rejection.

**Inference:** Distinct simulation seeds do not replace original-donor partitioning.

**Recommendation:** Use manifest-specific intersection; report coverage with conditional accuracy.

## L58 - M-PACT: deep-learning classification of pediatric CNS tumors from cell-free DNA methylation profiles.

Source: https://datacatalog.ccdi.cancer.gov/resource/M-PACT;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L59 - Artificial intelligence-generated synthetic data for cancer research and clinical trials.

Source: ../papers/s41568-026-00912-4.pdf; Sections on molecular generators, validation, and open challenges throughout

Status: Newly available 2026-09-10; methods-relevant sections reviewed

Delta category: CONFIRMS CURRENT PLAN

**Documented fact:** The review treats fidelity, utility, and privacy as separate dimensions and notes limited standardization.

**Inference:** A three-day generator adds several new claims without resolving cross-cancer transfer.

**Recommendation:** Keep augmentation gated, with untouched real validation only.

## L60 - In silico generation of synthetic cancer genomes using generative AI.

Source: ../papers/2025_In_silico_generation_of_synthetic_cancer_genomes_using_generative_AI.pdf; pp. 1-4,12

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Factorized generator ensemble evaluates genomic patterns and downstream classification; subclonal representation limited.

**Inference:** Joint CpG-intensity-allelic-CN covariance is not established.

**Recommendation:** Post-hackathon engineering inspiration only.

## L61 - Multi-omics data integration and drug screening of AML cancer using Generative Adversarial Network.

Source: ../papers/1-s2.0-S1046202324001075-main.pdf; Methods pp. 4-8; Results pp. 8-12

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** omicsGAN was evaluated on 173 TCGA AML cases with internal splits and repeated trials.

**Inference:** An internal AML AUC gain does not demonstrate pediatric or cancer-held-out HRDsum transfer.

**Recommendation:** Defer GAN/VAE/diffusion; any later augmentation must use identical untouched real validation.

## L62 - Performance analysis of data resampling on class imbalance and classification techniques on multi-omics data for cancer classification.

Source: ../papers/2024_Performance_analysis_of_data_resampling_on_class_imbalance_and_classification_te.pdf; pp. 1-5

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** SMOTE benefits reported on tumor-vs-normal tasks; reviewed text does not establish complete train-fold nesting.

**Inference:** Endpoint ease and unresolved leakage preclude borrowing accuracy claims.

**Recommendation:** Real-only primary; nested augmentation only.

## L63 - Opportunities and Challenges of Synthetic Data Generation in Oncology.

Source: ../papers/jacobs-et-al-2023-opportunities-and-challenges-of-synthetic-data-generation-in-oncology.pdf; Review sections on synthetic-data approaches, opportunities, and validation challenges

Status: Newly available 2026-09-10; relevance review completed

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** The article reviews opportunities and challenges for synthetic data generation in oncology.

**Inference:** Complex representation learning adds tuning and leakage risk without evidence of better cancer-held-out HRD prediction.

**Recommendation:** Defer neural representation models until the baseline and transfer benchmark are fixed.

## L64 - Genome-wide methylome modeling via generative AI incorporating long- and short-range interactions.

Source: ../papers/2025_Genome_wide_methylome_modeling_via_generative_AI_incorporating_long_and_short_ra.pdf; pp. 1-2,7-9

Status: Local abstract and relevant methods/results reviewed; supplementary gaps retained

Delta category: Not established in reviewed material

**Documented fact:** Imputes methylation using local/long-range context; millions of windows are not patients.

**Inference:** Imputation quality is not evidence of HRD transfer or label fidelity.

**Recommendation:** Defer; train-fold median imputation and dropout tests for MVP.

## L65 - Deep learning detects genetic alterations in cancer histology generated by adversarial networks.

Source: ../papers/The Journal of Pathology - 2021 - Krause - Deep learning detects genetic alterations in cancer histology generated by.pdf; Methods pp. 3-7; Results pp. 7-11

Status: Newly available 2026-09-10; full relevant methods/results reviewed

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** Combined real and synthetic histology training modestly improved held-out real-image AUC.

**Inference:** The evaluation design transfers to KIDS26; the generator architecture does not.

**Recommendation:** Any later augmentation comparison must use identical untouched real validation sets.

## L66 - Synthesis of diagnostic quality cancer pathology images by generative adversarial networks.

Source: https://pubmed.ncbi.nlm.nih.gov/32686118/;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L67 - Integrated genomic analyses of ovarian carcinoma.

Source: https://pubmed.ncbi.nlm.nih.gov/?term=Integrated+genomic+analyses+of+ovarian+carcinoma;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L68 - Comprehensive molecular portraits of human breast tumours.

Source: https://pubmed.ncbi.nlm.nih.gov/?term=Comprehensive+molecular+portraits+of+human+breast+tumours;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L69 - The Cancer Genome Atlas Pan-Cancer analysis project.

Source: https://pubmed.ncbi.nlm.nih.gov/?term=The+Cancer+Genome+Atlas+Pan-Cancer+analysis+project;

Status: Unavailable locally; metadata only; no methods inferred

Delta category: Not established in reviewed material

**Documented fact:** Listed in supplied map; full text unavailable locally.

**Inference:** No methodological inference from title.

**Recommendation:** Verify relevant methods before relying on this reference.

## L70 - Genomic signatures of homologous recombination deficiency predict response to PARP inhibition in chordoma

Source: ../papers/41467_2019_Article_9633.pdf; Main text pp. 1-10; Methods pp. 10-15

Status: Newly available 2026-09-10; locally present but absent from the original 69-row matrix

Delta category: INTERESTING BUT NOT HACKATHON-CRITICAL

**Documented fact:** Eleven advanced chordomas had matched tumor-normal sequencing; LOH regions exceeded 15 Mb and LST counted breaks between segments exceeding 10 Mb after small-segment filtering.

**Inference:** Rare-tumor relevance does not validate methylation-derived canonical HRDsum.

**Recommendation:** Use as biological motivation, not as evidence to change the primary model.