#!/bin/bash

#SBATCH --job-name=sims_binary_hist
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=600m
#SBATCH	--array=[1876-2500]
#SBATCH --output=/pine/scr/a/n/aniland/strapp_survival/Sims/Log/%a.out
#SBATCH --error=/pine/scr/a/n/aniland/strapp_survival/Sims/Error/%a.err

## add R module
module add gcc/11.2.0
module add r/4.1.0

R CMD BATCH --no-restore /proj/ibrahimlab/strapp_survival/Sims/R/binary_hist_sims.R /pine/scr/a/n/aniland/strapp_survival/Sims/Rout/array_$SLURM_ARRAY_TASK_ID.Rout