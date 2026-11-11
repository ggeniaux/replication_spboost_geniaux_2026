#!/bin/bash -l
# mstop-CV sub-experiment (see DESIGN.md): one task = one replication
# (= 4 cell x rho0 pairs), single-core.
#   cd /path/to/replication_scripts/e3prime && sbatch run_mstopcv.sbatch
#SBATCH -J e3mstopcv
#SBATCH --array=1-250%100
#   (site-specific options: partition, account, excluded nodes -> #SBATCH lines here)
#SBATCH -o slurm_out/%x_%a.out
#SBATCH -e slurm_out/%x_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=4G

cd "$SLURM_SUBMIT_DIR"
mkdir -p out_mstopcv slurm_out
module load ${R_MODULE:-R/4.5.1}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(hostname) rep=$SLURM_ARRAY_TASK_ID $(date)"
REPS="$SLURM_ARRAY_TASK_ID" NCORES=1 OUT=out_mstopcv Rscript run_mstopcv.R
echo "DONE rep=$SLURM_ARRAY_TASK_ID $(date)"
