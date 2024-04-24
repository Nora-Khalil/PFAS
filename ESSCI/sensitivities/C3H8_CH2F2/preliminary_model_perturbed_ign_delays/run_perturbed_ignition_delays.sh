#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=5:00:00
#SBATCH --job-name=perturbed_final_ignition_delays
#SBATCH --error=error_%x_%a.rmg.slurm.log
#SBATCH --output=output_%x_%a.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=short
#SBATCH --array=1-3

#values=("1818.18_0" "1818.18_4" "1467.18_0" "1467.18_4" "1292.51_0" "1292.51_4")
values=("1818.18_0.1" "1467.18_0.1" "1292.51_0.1")


index=$SLURM_ARRAY_TASK_ID-1

tuple="${values[$index]}" 

python perturbed_ignition_delays.py $tuple
