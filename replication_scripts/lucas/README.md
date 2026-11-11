# Lucas County workstream (empirical section of paper 1)

Scope (decisions 2026-11-05): estimation + diagnostics only,
no out-of-sample prediction. A priori menu of W {4-NN, 10-NN,
distance band}. Showcase: the "vintage" curve of `yrbuilt`.

## Execution order

```
Rscript 01_prep_lucas.R                      # lucas_prep.rds (data + W menu)
SMOKE=4000 Rscript 02_fit_lucas.R            # pipeline test (~min)
Rscript 02_fit_lucas.R                       # full run n=25357
Rscript 03_tables_figs_lucas.R               # paper tables + figures
```

## 02_fit_lucas.R — environment variables

- `WSEL=knn4,knn10,dband`: subset of the W menu.
- `MSTOP_MAX=2000`, `KSPEC=100`, `MU_GRID=0.01,1,100,1000,10000`.
- `SMOKE=<n>`: subsample for pipeline test (0 = full).
- `SKIP_BSPA=1`: skips the BSPA arm (the longest, mstop CV).
- `OUT=out_lucas`: RDS folder (`lucas_<W>.rds`).

## Contents of each `lucas_<W>.rds`

`free` / `free_cf` (free surface ± crossfit), `mu_path` (μ_s path),
`mu_reml` (REML indifference), `mu_star` + `star_cf` (plateau + crossfit,
final specification with var_rho / var_hc0 / var_shac), `bspa`
(selection, spatial-block CV mstop), `pik_free` / `pik_star` (π_k, 25
modes, shared eigenvectors).

## Points of attention

- `dband`: radius = smallest distance leaving zero isolates; check
  the max degree printed by 01 (risk of a dense W if a rural sale is
  very isolated). If max degree ≫ 100, discuss a cap before
  freezing the section.
- BSPA n=25k with mstop CV capped at 2000: expect several hours per W
  sequentially; run `SKIP_BSPA=1` first, then BSPA per W (e.g.
  on the HPC cluster, one job per W).
- `03_tables_figs_lucas.R` to be finalized once the RDS are complete:
  tab_lucas_rho.tex (ρ̂ × W × spec, 3 standard errors), fig_lucas_mu_path,
  fig_lucas_vintage (yrbuilt curve of the final spec), fig_lucas_pik.

## Map basemap

`lucas_basemap.rds` (county outline + water bodies projected in km in the
data CRS) is produced from the US Census TIGER files
(cb_2023_us_county_500k and tl_2023_39095_areawater, public domain,
downloaded into `tiger_raw/`, not versioned) by `00_make_basemap.R`.
