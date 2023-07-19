#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --job-name=PFPA
#SBATCH --error=rmg.slurm.lo_mg
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=west

source activate rmg_env
python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py
