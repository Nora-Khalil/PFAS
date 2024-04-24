#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=10-00:00:00
#SBATCH --job-name=2nd_try_autogenerate_H
#SBATCH --error=error_longer_run.slurm.log
#SBATCH --output=output_longer_run.slurm.log
#SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
##SBATCH --array=1
#SBATCH --partition=west


python BM_tree_fitting.py