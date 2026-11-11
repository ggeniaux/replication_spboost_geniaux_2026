# Replication scripts

Replication material for "Flexible nonlinear spatial autoregressive
models under spatial confounding: a gradient boosting approach with
closed-form estimation". All estimators are implemented in the
`spboost` R package (version 0.7.0), which these scripts assume is
installed; its source package is in `../package/`.

## Contents

### `e3prime/` — Monte Carlo designs (Section 3, Appendices)

- `DESIGN.md`: frozen factorial design of experiment E3' (axes S, C,
  F, rho, W) and of the E1/E2 campaigns.
- `run_e1e2.R`, `run_e3prime.R`, `run_mstopcv.R`: one-replication
  runners (one task = one design cell x one replication; seeds are set
  per cell for exact reproducibility). The `.sbatch`/`.sge` files are
  Slurm/SGE array-job wrappers for an HPC cluster; `reproduce.sh`
  chains the full campaign.
- `aggregate_*.R`: aggregation of raw per-replication outputs into the
  cell-level summaries `*_cells.csv` (included here, so that all
  tables and figures can be regenerated without re-running the
  campaigns).
- `gen_tab_e1e2.R`, `gen_tab_e3p.R`: LaTeX tables of the paper
  (Tables E1/E2, Table 7).
- `fig_e3prime_results.R`: Figures 2-3 (bias heatmaps, size barplot).
- `mc_sandwich.R`, `mini_crossfit.R`, `mini_kspec.R`: shared helpers
  (HC0/SHAC sandwich, spatial cross-fitting, spectral diagnostic).

### `lucas/` — Lucas County application (Section 4)

Numbered pipeline:

- `00_make_basemap.R`: county basemap from public US Census TIGER
  files (downloaded at run time).
- `01_prep_lucas.R`: data preparation from the `house` dataset of the
  `spData` R package (public).
- `02_fit_lucas.R`, `02b_fit_knn4_reml.R`: estimation under the three
  W matrices (4-NN, 10-NN, distance band), spectral-penalty path,
  cross-fitting, REML variant.
- `03_tables_figs_lucas.R`: Tables 8-9 and Figures 4-7.
- `04_blup_holdout.R` (+ `run_blup.sbatch`): hold-out prediction
  exercise.
- `05_extras_lucas.R`: impact multipliers (Hutchinson trace
  estimator) and auxiliary quantities.
- `out_lucas/`, `out_lucas_ext/`, `out_blup/`: fitted-model summaries
  (small `.rds`), included so tables and figures can be regenerated
  directly via `03_tables_figs_lucas.R` and `05_extras_lucas.R`
  without refitting.

### `fig_pi_k_E1_E2.R` — spectral diagnostic (Appendix B)

Self-contained script (simulates its own data with `spboost::dgp()`):
Figure B.10 (`pi_k` for designs E1, E2 and E1 with a geographic smoother),
plus the `spectral_overlap` screening and the dictionary approximation
error printed to the console. Output in `$OUT_DIR` (default `output/`).

## Reproducing from scratch

1. Install `spboost` 0.7.0 and its dependencies (R >= 4.2).
2. Lucas County: run `lucas/00` to `lucas/05` in order (a few hours on
   a laptop; `04` benefits from a cluster).
3. Monte Carlo: `e3prime/reproduce.sh` runs the full E1/E2/E3'
   campaigns (approximately 60-80 core-hours; designed for a
   Slurm/SGE array). The included `*_cells.csv` summaries allow
   regenerating every table and figure without this step.

Data: only simulated data and public data (spData `house` dataset, US
Census TIGER files) are used. No credentials or private resources are
required.
