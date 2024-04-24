#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=t_d_fluorine_lib
#SBATCH --error=error.slurm.log
#SBATCH --output=output.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-4
#SBATCH --partition=west


flame_scripts=(flame_speed_calc_range_1.py flame_speed_calc_range_2.py flame_speed_calc_range_3.py flame_speed_calc_range_4.py)

index=$SLURM_ARRAY_TASK_ID-1

file_name="${flame_scripts[$index]}"  


source activate cantera_env
python "$file_name"

