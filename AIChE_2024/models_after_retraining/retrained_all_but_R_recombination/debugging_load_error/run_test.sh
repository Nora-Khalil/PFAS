#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=R_Recombination
#SBATCH --error=R_Recombination.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=short

python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input_test.py
