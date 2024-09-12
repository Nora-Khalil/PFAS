#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=noFilt_PFPA_lowtol
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --partition=short

python-jl /home/khalil.nor/Code/RMG-Py/rmg.py input.py

