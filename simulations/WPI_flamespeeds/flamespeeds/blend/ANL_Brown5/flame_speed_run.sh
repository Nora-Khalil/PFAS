#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=A_B5_blend
#SBATCH --error=error_fc.%a.slurm.log
#SBATCH --output=output_fc.%a.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-10
#SBATCH --partition=short

index=$SLURM_ARRAY_TASK_ID-1

python flamespeeds_slurm-array.py $index

