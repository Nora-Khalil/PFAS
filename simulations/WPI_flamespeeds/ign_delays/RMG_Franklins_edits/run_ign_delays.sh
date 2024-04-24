#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --job-name=automated_ign_delays
#SBATCH --error=error_fc_%a.slurm.log
#SBATCH --output=output_fc_%a.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-80
#SBATCH --partition=express

index=$SLURM_ARRAY_TASK_ID-1

python sbatch_ign_delays.py $index