# Team and Roles

- **Team name:** KIDS26 Team 12 — Methyl-HRD
- **Team leads:** [Evan Savage (@esavage111)](https://github.com/esavage111) and [Susanna Downing (@sdowning12)](https://github.com/sdowning12)
- **Communication channel:** [Team 12 Slack](https://stjudebiohackathon.slack.com/archives/C0BSC28M3U6)
- **General channel:** [St. Jude Biohackathon general Slack](https://stjudebiohackathon.slack.com/archives/C04JD4M3TCM)
- **Mentor/support contact:** Evan Savage (@esavage111)
- **Project question/problem:** Can adult TCGA methylation predict allele-aware reference HRDsum in unseen cancer types and transfer without refitting to eligible pediatric high-grade glioma samples?
- **Expected output:** Auditable data and feature artifacts, a leakage-safe cancer-held-out elastic-net result with null/OOD reporting, and a precomputed Shiny demonstration.
- **Tools and stack:** Python 3.10+, R/`glmnet`/`renv`, optional Bioconductor/conumee2, Shiny, GitHub, and Slack.

## Roles

| Person | Role | Main responsibility | Backup or support needed |
|---|---|---|---|
| Evan Savage (@esavage111) | Co-lead; data and integration lead; mentor/support contact | Repository coordination, TCGA acquisition/provenance, cohort construction, pipeline integration, reproducibility, and blocker resolution | Susanna reviews scientific gates and claims |
| Susanna Downing (@sdowning12) | Co-lead; methods and validation lead | Endpoint/partition freeze, leakage review, PBTP readiness coordination, interpretation, figures, and presentation | Evan supports implementation and data provenance |

Roles overlap. At the first check-in, the leads should assign named owners for model execution, independent result review, CNV gate, Shiny/demo, and presentation based on the participants present. No unlisted person or role is invented in this document.
