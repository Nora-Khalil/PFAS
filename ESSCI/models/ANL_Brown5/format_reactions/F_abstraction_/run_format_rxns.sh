#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=format_F_rxns
#SBATCH --error=error.slurm.log
#SBATCH --output=output.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


python format_F_abstraction_reactions.py