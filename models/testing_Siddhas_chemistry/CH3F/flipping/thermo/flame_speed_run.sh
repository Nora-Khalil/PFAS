#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=siddhas_flipped_thermo
##SBATCH --error=./errors/error_fc_%a.slurm.log
##SBATCH --output=./outputs/output_fc_%a.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-17
#SBATCH --partition=short

list_of_ctis=(copy_12.cti copy_128.cti copy_161.cti copy_19.cti copy_3313.cti copy_37.cti copy_4.cti copy_42.cti copy_43.cti copy_44.cti copy_45.cti copy_67.cti copy_68.cti copy_7.cti copy_75.cti copy_76.cti copy_8.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_ctis[$index]}"  


source activate cantera_env
python flame_speed_calc.py $folder_name
