#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=16-00:00:00
#SBATCH --job-name=PFOS_H2O
#SBATCH --error=rmg.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --partition=west

source activate rmg_env
python-jl /home/khalil.nor/Code/RMG-Py/rmg.py -n 5 -p input.py

