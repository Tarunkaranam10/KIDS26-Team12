# KIDS26 Biohackathon Project Plan

This repository contains Team 12's three-day biohackathon project: a reproducible test of whether DNA methylation can predict independently measured homologous recombination deficiency genomic-scar burden across cancer types and transfer to pediatric high-grade glioma.

> **Team leads:** Start with the [team lead checklist](project-management/CHECKLIST.md) before the event or during the first team meeting.

## Project Profile

- **Project name:** KIDS26 Methyl-HRD: Cross-cancer prediction of genomic-scar burden
- **Question, problem, or opportunity:** Can a leakage-safe model trained on adult TCGA 450K methylation predict continuous, independently measured reference HRDsum in cancer types it has not seen and then transfer, without refitting, to eligible PBTP pediatric high-grade glioma samples?
- **Data, inputs, or evidence:** Public TCGA/PanCanAtlas 450K methylation, PanImmune HRD-LOH/LST/TAI/HRDsum labels, subtype/QC/purity/ploidy metadata, versioned HM450/EPIC-v1 annotations, 60 locally reviewed papers, and approved de-identified PBTP methylation plus allele-aware genomic reference labels if available.
- **Expected output:** An auditable TCGA cohort, frozen 450K/EPIC feature bridge, nested cancer-aware elastic-net baseline, held-out-cancer and null-comparator results, uncertainty/OOD reporting, and a Shiny demonstration using public, synthetic, or approved de-identified precomputed results.
- **Tools and stack:** R 4.5.0 or Python 3.10+ for acquisition, matching, QC, and reproducibility; R with `glmnet`, `data.table`, `jsonlite`, and `renv` for modeling; Bioconductor/minfi and conumee2 only for the gated intensity/CNV branch; Shiny for presentation; Git/GitHub for collaboration; Slack for team communication.
- **Team leads:** [Evan Savage (@esavage111)](https://github.com/esavage111) and [Susanna Downing (@sdowning12)](https://github.com/sdowning12)
- **Team members and roles:** [Team and roles](project-management/team.md)
- **Communication:** [Team 12 Slack channel](https://stjudebiohackathon.slack.com/archives/C0BSC28M3U6); [Biohackathon general channel](https://stjudebiohackathon.slack.com/archives/C04JD4M3TCM)

The stack and roles may be narrowed during the event. Optional methods cannot displace the auditable real-data baseline and unseen-cancer validation.

## Vision and Mission

- **Vision:** Enable a scientifically defensible, single-sample methylation assay to estimate genomic-scar burden in pediatric tumors while exposing uncertainty, domain shift, and the limits of methylation-derived copy-number information.
- **Mission:** During the biohackathon, construct and freeze the public TCGA training cohort and feature bridge, fit a leakage-safe elastic-net model, test cross-cancer generalization against null comparators, and demonstrate approved precomputed results without making unsupported clinical or treatment claims.

## About

Homologous recombination deficiency can leave persistent genomic scars associated with DNA-repair defects and treatment sensitivity. Canonical HRDsum combines HRD-LOH, large-scale state transitions, and telomeric allelic imbalance measured from allele-aware genomic data. Methylation arrays provide rich CpG profiles and can support total-copy-number inference from raw intensities, but ordinary total CN cannot identify copy-neutral LOH or prove telomeric allelic imbalance.

Team 12 therefore treats methylation as a predictor of an independently measured continuous reference HRDsum. The primary experiment uses shared HM450/EPIC CpGs and elastic-net regression. Cancer-type-held-out validation tests whether performance reflects transferable HRD biology rather than tissue identity. Total-CN features and fusion are secondary gated comparisons; synthetic augmentation remains a stretch experiment after the real-data result is stable.

This is a research prototype. It is not a clinical HRD assay, functional repair test, or treatment recommendation.

## Roadmap and Milestones

| When | Focus | Expected outcome |
|---|---|---|
| Day 1 | Verify downloaded inputs; freeze the sample audit, 450K/EPIC bridge, roles, and analysis protocol; run R smoke tests and the first real elastic-net/null comparison | An auditable eligible cohort and at least one reproducible cancer-held-out result by early Day 2 |
| Day 2 | Complete development leave-one-cancer-type-out evaluation, select and freeze the model, then apply it to locked CNS and eligible PBTP data; gate CNV/fusion | A working transfer result or clear evidence that the hypothesized transfer does not hold |
| Day 3 | Finish uncertainty/OOD and interpretation by noon; stabilize figures, Shiny demo, documentation, presentation, and rehearsal | A reproducible demonstration or negative result with methods, limitations, provenance, and next steps |

The goal is a clear, honest, useful result that the team can explain and others can build on. Optional modeling stops before it threatens the primary result or the second half of Day 3.

## Start Here

- [Detailed technical README](README.md)
- [Project plan](project-management/project-plan.md)
- [Three-day runbook](docs/05_THREE_DAY_HACKATHON_RUNBOOK.md)
- [Pre-hackathon checklist](docs/04_PRE_HACKATHON_CHECKLIST.md)
- [Data acquisition and cohort construction](docs/03_DATA_ACQUISITION.md)
- [Verification status](docs/14_VERIFICATION_STATUS.md)
- [Literature reconciliation](docs/15_NEW_LITERATURE_RECONCILIATION.md)
- [PBTP readiness specification](docs/18_PBTP_READINESS.md)
