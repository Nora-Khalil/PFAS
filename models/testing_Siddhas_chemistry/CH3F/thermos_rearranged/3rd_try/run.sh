#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=6-00:00:00
#SBATCH --job-name=CH3F_thermos_rearranged_family
#SBATCH --error=error.rmg.slurm.log
#SBATCH --output=output.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=west

source activate rmg_env
python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py