# Replication material — spboost (Geniaux, 2026)

Replication material for

> Geniaux, G. (2026). *Flexible nonlinear spatial autoregressive models
> under spatial confounding: a gradient boosting approach with closed-form
> estimation.* Preprint.

All estimators are implemented in the `spboost` R package (version 0.7.0),
whose source package is included here.

## Contents

| Path | Content |
|---|---|
| `package/spboost_0.7.0.tar.gz` | source package of `spboost` 0.7.0 (GPL (>= 2)) |
| `replication_scripts/e3prime/` | Monte Carlo designs E1/E2/E3' (Section 3 and appendices) |
| `replication_scripts/lucas/` | Lucas County, Ohio, housing application (Section 4) |
| `replication_scripts/README.md` | detailed description of every script |

The replication-level aggregates `*_all.rds` (up to 33 MB) are distributed
only in this repository, not in the ancillary files of the arXiv/HAL
preprint.

Small precomputed summaries (`*_cells.csv`, fitted-model `.rds` files) are
included so that every table and figure of the paper can be regenerated
without re-running the Monte Carlo campaigns or refitting the models.

## Quick start

Requirements: R (>= 4.2) and a C++ toolchain.

```sh
R CMD INSTALL package/spboost_0.7.0.tar.gz
```

Install its CRAN dependencies first, plus those used by the scripts:

```r
install.packages(c("Matrix", "mboost", "mgcv", "mgwrsar", "Rcpp",
  "RcppEigen", "sf", "MASS", "data.table", "xgboost", "caret",
  "doParallel", "foreach", "nabor", "earth",
  "RSpectra", "spdep", "spatialreg", "spData"))
```

- **Lucas County** — run `replication_scripts/lucas/00_*.R` to `05_*.R` in
  order (a few hours on a laptop; `04` benefits from a cluster). Tables and
  figures alone: `03_tables_figs_lucas.R` and `05_extras_lucas.R`, which
  read the included `.rds` summaries.
- **Monte Carlo** — `replication_scripts/e3prime/reproduce.sh [smoke|pilot|full|minis]`
  runs the campaigns locally (`full` ≈ 60–80 core-hours; Slurm/SGE array
  wrappers `*.sbatch` / `*.sge` are provided). Tables and figures alone:
  `gen_tab_e1e2.R`, `gen_tab_e3p.R` and
  `PREFIX=e3prime_full Rscript fig_e3prime_results.R`, from the included
  cell summaries (`*_cells.csv`) and replication-level aggregates
  (`e1e2_all.rds`, `e3prime_full_all.rds`, `e3prime_mstopcv_all.rds`, as
  produced by the `aggregate_*.R` scripts).

Tables and figures are written to `output/` (override with the `OUT_DIR`
environment variable). Seeds are set per design cell and replication for
exact reproducibility.

## Data

- House sales: dataset `house` of the public R package `spData`.
- County boundary and water features: public US Census TIGER files,
  downloaded at run time by `lucas/00_make_basemap.R` (the only step that
  needs network access).

## Licence

Code and package: GNU General Public License, version 2 or (at your option)
any later version — see `LICENSE`.

## Citation

See `CITATION.cff`. Please cite the paper above when using this material.
