#!/bin/bash

#SBATCH --job-name=sims_binary_hist
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=600m
#SBATCH	--array=[1-10]
#SBATCH --output=/proj/ibrahimlab/strapp_survival/Sims/Results/Log/slurmLogFiles%a.out
#SBATCH --error=/proj/ibrahimlab/strapp_survival/Sims/Results/Error/%a.err

## add R module
module add gcc/11.2.0
module add r/4.1.0

R CMD BATCH --no-restore /proj/ibrahimlab/strapp_survival/Sims/R/binary_hist_sims.R /proj/ibrahimlab/strapp_survival/Sims/Results/Rout/array_$SLURM_ARRAY_TASK_ID.Rout