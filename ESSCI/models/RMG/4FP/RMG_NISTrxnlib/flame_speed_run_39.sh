#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=7:00:00
#SBATCH --job-name=fs_39_NISTrxnlib
#SBATCH --error=error_fc_39.slurm.log
#SBATCH --output=output_fc_39.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


python flame_speed_calc_39.py

