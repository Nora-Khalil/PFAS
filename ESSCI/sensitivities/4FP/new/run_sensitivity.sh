#!/bin/bash
#SBATCH --job-name=s_RMG

#SBATCH --partition=short
#SBATCH --time=15:00:00
#SBATCH --exclude=c5003
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=4



python sensitivity.py 
