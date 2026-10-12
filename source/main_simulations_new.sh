#!/bin/bash
#SBATCH --array=1-81%81
#SBATCH --job-name=main_simulations_new
#SBATCH --partition=wrobel
#SBATCH --output=main_simulations_new_%a.out
#SBATCH --error=main_simulations_new_%a.err

# 81 scenarios = 3 kappa x 3 N x 3 tau x 3 censoring rates.
# The R script refuses any index outside 1..nrow(params), so a grid/array
# mismatch fails loudly instead of silently skipping scenarios.

module purge
module load R/4.4.0

JOBID=$SLURM_ARRAY_TASK_ID
Rscript main_simulations_new.R $JOBID
