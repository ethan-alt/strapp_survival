#!/bin/bash

#SBATCH --job-name=strapp_analysis
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250m
#SBATCH	--array=[1-65]
#SBATCH --output=/proj/ibrahimlab/strapp_survival/Analysis/Results/Log/slurmLogFiles%a.out
#SBATCH --error=/proj/ibrahimlab/strapp_survival/Analysis/Results/Error/%a.err

## add R module
module add gcc/11.2.0
module add r/4.1.0

R CMD BATCH --no-restore /proj/ibrahimlab/strapp_survival/Analysis/R/01_data_analysis.R /proj/ibrahimlab/strapp_survival/Analysis/Results/Rout/analysis_$SLURM_ARRAY_TASK_ID.Rout