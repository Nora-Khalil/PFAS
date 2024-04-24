#!/bin/sh

#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=fc_CH2F2
#SBATCH --error=error_fc.slurm.log
#SBATCH --output=output_fc.slurm.log
##SBATCH --cpus-per-task=5
##SBATCH --mem-per-cpu=8Gb
##SBATCH --ntasks=1 
#SBATCH --array=1-60
#SBATCH --partition=short


list_of_ctis=(copy_184_chem.cti copy_184_chem0006.cti copy_184_chem0007.cti copy_184_chem0008.cti copy_184_chem0009.cti copy_184_chem0010.cti copy_184_chem0011.cti copy_184_chem0012.cti copy_184_chem0013.cti copy_184_chem0014.cti copy_184_chem0015.cti copy_184_chem0016.cti copy_184_chem0017.cti copy_184_chem0018.cti copy_184_chem0019.cti copy_184_chem0020.cti copy_184_chem0021.cti copy_184_chem0022.cti copy_184_chem0023.cti copy_184_chem0024.cti copy_184_chem0025.cti copy_184_chem0026.cti copy_184_chem0027.cti copy_184_chem0028.cti copy_184_chem0029.cti copy_184_chem0030.cti copy_184_chem0031.cti copy_184_chem0032.cti copy_184_chem0033.cti copy_184_chem0034.cti copy_184_chem0035.cti copy_184_chem0036.cti copy_184_chem0037.cti copy_184_chem0038.cti copy_184_chem0039.cti copy_184_chem0040.cti copy_184_chem0041.cti copy_184_chem0042.cti copy_184_chem0043.cti copy_184_chem0044.cti copy_184_chem0045.cti copy_184_chem0046.cti copy_184_chem0047.cti copy_184_chem0048.cti copy_184_chem0049.cti copy_184_chem0050.cti copy_184_chem0051.cti copy_184_chem0052.cti copy_184_chem0053.cti copy_184_chem0054.cti copy_184_chem0055.cti copy_184_chem0056.cti copy_184_chem0057.cti copy_184_chem0058.cti copy_184_chem0059.cti copy_184_chem0060.cti copy_184_chem0061.cti copy_184_chem0062.cti copy_184_chem0063.cti copy_184_chem_annotated.cti)

index=$SLURM_ARRAY_TASK_ID-1

folder_name="${list_of_ctis[$index]}" 


source activate cantera_env
python flame_speed_calc-one-point.py $folder_name


