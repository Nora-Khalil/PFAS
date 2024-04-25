#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=25-00:00:00
#SBATCH --job-name=4FP_retrained
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=west

python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py

