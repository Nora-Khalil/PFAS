#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=15-00:00:00
#SBATCH --job-name=A_B5_pure_propane
#SBATCH --error=error_fc_loop.slurm.log
#SBATCH --output=output_fc_loop.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1-40
#SBATCH --partition=west

python flamespeeds.py