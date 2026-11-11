#!/bin/bash
# Reproduction of the E3' campaign WITHOUT a cluster (see DESIGN.md).
# Usage :  bash reproduce.sh [smoke|pilot|full|minis]
#   smoke : 5 replications   (~15 min on 8 cores)  - pipeline check
#   pilot : 100 replications (~3 h on 8 cores)
#   full  : 1000 replications (paper campaign; ~30 h on 8 cores,
#           prefer the cluster: run_e3prime.sbatch / run_e3prime.sge)
#   minis : the two mini-experiments that fixed the design
#           (mini_kspec.R: kspec and per-term spx; mini_crossfit.R: cross-fitted
#           CFE), appendix tables.
# Variables: NCORES (default: cores-2), OUT (default: out_<mode>).
# Prerequisites: the spboost package from THIS commit installed
#   (R CMD INSTALL . at the repository root), RSpectra, mboost, mgcv.
# Partial restart: already computed pairs (RDS present) are skipped.
set -euo pipefail
cd "$(dirname "$0")"
MODE=${1:-smoke}
NCORES=${NCORES:-$(Rscript -e 'cat(max(1L, parallel::detectCores() - 2L))')}

Rscript -e 'stopifnot(requireNamespace("spboost", quietly = TRUE));
  stopifnot("crossfit" %in% names(formals(spboost::GAM_SAR_CFE)));
  cat("spboost", as.character(packageVersion("spboost")),
      "(with crossfit) from", find.package("spboost"), "\n")'

case "$MODE" in
  smoke) REPS="1:5" ;;
  pilot) REPS="1:100" ;;
  full)  REPS="1:1000" ;;
  minis)
    Rscript mini_kspec.R
    Rscript mini_crossfit.R
    exit 0 ;;
  *) echo "unknown mode: $MODE (smoke|pilot|full|minis)"; exit 1 ;;
esac

OUT=${OUT:-out_$MODE}
echo "== E3' $MODE : REPS=$REPS, $NCORES cores, outputs $OUT =="
REPS=$REPS NCORES=$NCORES OUT=$OUT RHOS="0,0.2,0.6" Rscript run_e3prime.R
OUT=$OUT Rscript aggregate_e3prime.R
echo "== done: e3prime_pilot_cells.csv / _traj.csv / _all.rds =="
