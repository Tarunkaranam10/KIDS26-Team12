# Reproducible environment status

Use Python >=3.10 for acquisition/testing (standard library only). `pypdf` is needed only to re-extract source PDFs, not to reproduce acquisition/model inputs. The live setup used bundled Python; record `python --version` on your analysis host.

Use R with `glmnet`, `data.table`, `jsonlite`, `shiny`, and `renv` for the baseline/demo. R/Bioconductor `sesame`, `minfi`, and conumee2 belong to the separate intensity branch. Pin the Bioconductor release compatible with the selected R and caller versions; do not invent a renv.lock without installing and testing.

```sh
Rscript scripts/setup.R
Rscript tests/smoke_model.R
Rscript -e "renv::snapshot(prompt=FALSE)"
```

Commit the resulting tested lockfile and `renv/activate.R` if generated; keep libraries ignored. `renv::restore()` is appropriate only after that lockfile exists. The supplied setup intentionally does not claim an already locked or tested R environment. No R runtime was found in this Windows task, so R execution is pending. Python tests and source compilation are recorded separately.

Run from repository root. Acquire/prepare data before fitting; the R script rejects engineering-only probe matrices and requires multi-cancer development data. Full historical matrices need streaming feature reduction; the R fitting script expects a bounded feature matrix, not the unfiltered 41.5 GB input. Seed 260910; save sessionInfo, model/input hashes and exact Git revision with every result.

A container is a post-validation packaging step. Do not use a container image tag as a substitute for a reproducible lockfile and successful smoke tests.
