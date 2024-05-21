#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=5-00:00:00
#SBATCH --job-name=R_recombination
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=west

python R_Recombination.py