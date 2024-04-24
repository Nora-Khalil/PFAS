#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=autogenerate_F
#SBATCH --error=error_F.slurm.log
#SBATCH --output=output_F.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=short


python BM_tree_fitting.py