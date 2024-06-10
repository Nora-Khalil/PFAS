#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=t_d_spec_core
#SBATCH --error=error_1.slurm.log
#SBATCH --output=output_1.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --partition=short 

python flame_speed_calc_range_1.py 

