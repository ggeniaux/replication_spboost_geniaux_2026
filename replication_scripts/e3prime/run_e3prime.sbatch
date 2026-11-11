#!/bin/bash -l
# (login shell: makes `module` available on the compute nodes)
# E3' campaign (spboost, the HPC cluster/slurm, frozen design DESIGN.md): one task =
# one replication (= 27 S x C x rho0 pairs, all arms, single-core).
#   cd /path/to/replication_scripts/e3prime && mkdir -p slurm_out && sbatch run_e3prime.sbatch
# Limit concurrent tasks (%N suffix of --array) according to site rules.
#SBATCH -J e3full
#SBATCH --array=1-1000%100
#   (site-specific options: partition, account, excluded nodes -> #SBATCH lines here)
#SBATCH -o slurm_out/%x_%a.out
#SBATCH -e slurm_out/%x_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=4G

cd "$SLURM_SUBMIT_DIR"
mkdir -p "${OUT:-out_full}" slurm_out
module load ${R_MODULE:-R/4.5.1}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(hostname) rep=$SLURM_ARRAY_TASK_ID $(date)"
REPS="$SLURM_ARRAY_TASK_ID" NCORES=1 OUT="${OUT:-out_full}" RHOS="${RHOS:-0,0.2,0.6}" N="${N:-1000}" CELLS="${CELLS:-}" Rscript run_e3prime.R
echo "DONE rep=$SLURM_ARRAY_TASK_ID $(date)"
