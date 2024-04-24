#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=perturbed_ignition_delays
#SBATCH --error=error_%x_%a.rmg.slurm.log
#SBATCH --output=output_%x_%a.slurm.log
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8Gb
#SBATCH --ntasks=1 
#SBATCH --array=1
#SBATCH --partition=short
#SBATCH --array=1-15

values=("1818.18_0" "1818.18_0.1" "1818.18_0.5" "1818.18_2" "1818.18_4" "1467.18_0" "1467.18_0.1" "1467.18_0.5" "1467.18_2" "1467.18_4" "1292.51_0" "1292.51_0.1" "1292.51_0.5" "1292.51_2" "1292.51_4")


index=$SLURM_ARRAY_TASK_ID-1

tuple="${values[$index]}" 


python perturbed_ign_delays.py $tuple
