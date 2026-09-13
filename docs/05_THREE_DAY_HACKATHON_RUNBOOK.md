# Three-day execution runbook

Times assume 09:00-17:00 with lunch 12:00-13:00. Roles can run independent work in parallel; this does not require extra software agents. No major modeling after Day 3 noon.

| Block | Owner | Input | Executable task | Output / success criterion | Fallback / stop-go |
|---|---|---|---|---|---|
| D1 09-10 | Data + reviewer | Downloads/manifests | Run build_master; audit label sums, groups, Ns and exclusions | Frozen eligible table + sample flow | No reliable labels: stop fitting, repair joins |
| D1 10-12 | Model | Eligible betas/allowlist | prepare_beta; R smoke tests; first grouped within-adult baseline | Reproducible fit/null metric table; all transforms inside folds | Downsize features/cancers; no global preprocessing shortcut |
| D1 10-12 | CNV | IDATs/normals/paired reference | Run pinned caller on 3-5 specimens; inspect signal blind to outcome | Caller and coordinate smoke test | No intensities/normals: cut CNV/fusion |
| D1 13-15 | Model + stats | Baseline dataset | Run development LOCO; capture tuning/fold IDs | First real held-out cancer predictions and nulls | Too few cancer groups: label within-cancer proof-of-concept |
| D1 13-15 | Data | Full matrix/intensity transfer | Finish staged download/QC; resolve replicate exclusions | Adequate per-cancer Ns and truth coverage | Freeze reduced cohort; document selection |
| D1 15-17 | Reviewer + all | Initial predictions/CNV | Inspect macro/per-cancer errors; permutation/null checks | Minimum viable scientific result or a documented failed hypothesis | Do not optimize on held-out CNS/PBTP |
| D2 09-11 | CNV + stats | Paired reference subset | Fixed-bin agreement, arms, breakpoints, purity strata | CNV benchmark with per-patient errors | Gate fails: report negative result; no CN predictors |
| D2 09-11 | Model | Non-CNS development | Finish nested tuning/LOCO; compare training-only mean | Stable primary architecture and tuning budget | Keep simple baseline; cut new models |
| D2 11-12 | Model | Paired CpG/CNV features | Fit B/C on same patients/splits | Incremental macro-MAE comparison | No paired gain: choose CpG-only |
| D2 13-14 | Stats/reviewer | Development comparisons | Select architecture, freeze preprocessing/features/model/calibration | Saved hashes, session, protocol and output schema | Any unresolved QC: defer external test |
| D2 14-16 | Independent evaluator | Locked CNS then permitted PBTP | Apply frozen model, no refit; patient-group metrics | Transfer evaluation including abstentions | No PBTP: CNS-only scope. Adaptation becomes separately labeled exploration |
| D2 16-17 | Model / demo | Stable results | Optional component analysis or one PCA-interpolation sensitivity; start precomputed Shiny | Optional comparison plus demo schema | Skip augmentation immediately if core unstable |
| D3 09-10 | Stats | Fixed predictions | Bootstrap patient/lineage metrics; interval coverage and OOD | Complete uncertainty/error tables | Report broad uncertainty; no unsupported probabilities |
| D3 10-11 | Figure/demo | Approved outputs | Validation scatter, per-cancer errors, sample flow, selected CN plot | Figures state N, assay and target | Synthetic fixtures clearly labeled if real output unavailable |
| D3 11-12 | Reviewer + lead | Code/results | Critical bug checks, rerun affected analysis only; freeze release | Final hashes and reproducibility instructions | Critical flaw: downgrade claims, preserve honest demo |
| D3 13-14 | Demo | Frozen approved outputs | Shiny display and failure-state rehearsal | Runs offline; no unapproved patient upload | Static plots/table fallback |
| D3 14-15 | Lead + science | Frozen conclusions | Slides: hypothesis, data, method, transfer, limits | One claim per result; distinguish negative findings | Cut optional panels |
| D3 15-16 | Team | Presentation | Rehearse timing, questions and handoff | Every claim traceable to a result/source | Fix presentation only |
| D3 16-17 | Team | Frozen release | Demonstration/presentation | Reproducible scientific story | No speculative model changes |

Cut order if late: diffusion/GAN/VAE -> RNA -> synthetic interpolation -> nonlinear trees -> components -> fusion -> CNV prediction. Preserve sample audit, simple baseline, unseen-cancer validation, honest uncertainty and reproducibility. If public data arrive late, a well-tested acquisition/matching pipeline and predeclared experiment are useful preparation but must not be called successful biological prediction.

## Hard milestone hierarchy

- **By early Day 2:** auditable cohort, frozen feature bridge, leakage-safe CpG elastic-net baseline, at least one held-out-cancer result, and a training-only-mean null comparator.
- **Strong result:** complete multi-cancer transfer audit, locked CNS evaluation, PBTP application only if authorized/ready, OOD/uncertainty reporting, and interpretation.
- **Stretch only:** CNV branch, fusion, component models, synthetic augmentation, and advanced Shiny features.

No new literature result changes the noon Day-3 experimentation stop. The laboratory-consensus evidence adds a mandatory check before slides: every displayed HRD value must name its source assay/caller and avoid a portable binary cutoff unless independently validated.
