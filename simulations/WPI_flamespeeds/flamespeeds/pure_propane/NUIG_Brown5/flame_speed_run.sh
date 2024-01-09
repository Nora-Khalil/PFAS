#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=N_B5_pure_propane
#SBATCH --error=error_fc_%a.slurm.log
#SBATCH --output=output_fc_%a.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-40
#SBATCH --partition=short

index=$SLURM_ARRAY_TASK_ID-1

python flamespeeds_sbatch_arr.py $index