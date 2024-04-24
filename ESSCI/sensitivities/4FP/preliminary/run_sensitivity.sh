#!/bin/bash
#SBATCH --job-name=s_pre_RMG

#SBATCH --partition=short
#SBATCH --time=5:00:00
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4



python sensitivity.py 
