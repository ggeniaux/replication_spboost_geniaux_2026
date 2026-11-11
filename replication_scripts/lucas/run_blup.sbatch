#!/bin/bash -l
#SBATCH --job-name=lucas_blup
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#   (site-specific options: partition, account, excluded nodes -> #SBATCH lines here)
#SBATCH --output=logs/blup_%a.out

# Lucas out-of-sample BLUP experiment (dispersed 10% test set, 4 BSPA arms).
# Prerequisites in /path/to/replication_scripts/lucas: 04_blup_holdout.R, lucas_prep.rds
# (light version: d + meta suffice, the W matrices are rebuilt by the script).
module load ${R_MODULE:-R/4.5.1}
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs out_blup
REPS="${SLURM_ARRAY_TASK_ID}" NREP=10 MSTOP_MAX=6000 OUT=out_blup \
  ARMS="${ARMS:-nos_rand,nos_spat,s_rand,s_spat}" Rscript 04_blup_holdout.R
