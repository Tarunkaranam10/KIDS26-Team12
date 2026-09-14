# Team and Roles

- **Team name:** Genomic Scars Pipeline: Assessment and Visualization Shiny (Team 12)
- **Team lead(s):** Evan Savage (@esavage111) and Susanna Downing (@sdowning12)
- **Communication channel:** [Add link]
- **Project question/problem:** Can a frozen methylation model predict independently measured HRDsum in a cancer type it has not seen, ultimately including pediatric high-grade glioma?
- **Expected output:** A methylation-based classifier that predicts genomic scar burden (HRDsum) trained on TCGA data and applied to a pediatric cancer set, delivered with a Shiny visualization application.
- **Tools and stack:** Methylation array data (CpG-level), elastic-net regression, R/Shiny for visualization, TCGA as training data, pediatric cancer cohort (incl. high-grade glioma) as the independent application set.

## Approach Summary

**Primary approach:** Technically harmonized autosomal CpGs → elastic-net regression → predicted reference HRDsum. Canonical HRDsum = HRD-LOH + LST + TAI requires allele-aware genomic information; methylation total-CN profiles cannot directly identify LOH/TAI. CNV-only and fusion predictors are gated secondary comparisons. The prototype is a research scar-burden predictor, not a validated functional-HRD or therapy-selection test.

## Roles

| Person | GitHub Handle | Role | Main responsibility | Backup or support needed |
| --- | --- | --- | --- | --- |
| Evan Savage | @esavage111 | Project Lead | Overall project direction, coordination, methylation modeling strategy | — |
| Susanna Downing | @sdowning12 | Co-Lead / Technical Contributor, Data Explorer, Domain Expert | Domain expertise on HRD biology, data exploration, elastic-net model development | Evan Savage |
| Tarun Venkat Sai Karanam | @Tarunkaranam10 | Technical Contributor, Data Explorer, Product Designer | Pipeline implementation, TCGA data exploration, Shiny app UI/UX | Sanika Naik |
| Sai Thulabandu | [confirm handle] | Data Explorer | TCGA/pediatric data exploration and QC | Talia Dalton |
| Kayode Raheem | @kayoderaheem | Technical Contributor, Data Explorer | Methylation preprocessing, model harmonization support | Ulofe Uduokhai |
| Sanika Naik | [confirm handle] | Product Designer, Data Explorer, Domain Expert | Shiny app design, HRD domain input, data exploration | Tarun Venkat Sai Karanam |
| Talia Dalton | @talidalton | Data Explorer | Pediatric cancer dataset exploration and QC | Sai Thulabandu |
| Ulofe Uduokhai | @c4usal | Technical Contributor, Data Explorer, Product Designer | Elastic-net model support, Shiny app development | Kayode Raheem |
| Eshwar P. Alicom | [confirm handle] | Domain Expert, Technical Contributor | HRD/genomic scar domain expertise, model validation input | Henry Dinh |
| Henry Dinh | @HenryDinh2005 | Domain Expert | Domain expertise on pediatric high-grade glioma and cancer genomics | Eshwar P. Alicom |

Roles can overlap. Revisit them when the project direction or stack changes.
