#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=t_d_fluorine_lib
#SBATCH --error=error_2.slurm.log
#SBATCH --output=output_2.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --partition=short


python flame_speed_calc_range_2.py 

