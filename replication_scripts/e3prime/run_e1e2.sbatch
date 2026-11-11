#!/bin/bash -l
# Rerun of E1/E2(i), restricted scope (DESIGN.md): one task = one
# replication (60 cells x estimators), single-core.
#   cd /path/to/replication_scripts/e3prime && sbatch run_e1e2.sbatch
#SBATCH -J e1e2
#SBATCH --array=1-1000%100
#   (site-specific options: partition, account, excluded nodes -> #SBATCH lines here)
#SBATCH -o slurm_out/%x_%a.out
#SBATCH -e slurm_out/%x_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=6G

cd "$SLURM_SUBMIT_DIR"
mkdir -p "${OUT:-out_e1e2}" slurm_out
module load ${R_MODULE:-R/4.5.1}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(hostname) rep=$SLURM_ARRAY_TASK_ID $(date)"
REPS="$SLURM_ARRAY_TASK_ID" NCORES=1 OUT="${OUT:-out_e1e2}" Rscript run_e1e2.R
echo "DONE rep=$SLURM_ARRAY_TASK_ID $(date)"
