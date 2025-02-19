#!/usr/bin/env bash
#SBATCH -A naiss2024-22-1647
#SBATCH -J mip
#SBATCH -n 64
#SBATCH -t 0-24:00:00
#SBATCH -p shared
#SBATCH --mem=80240



######################################################
module load bioinfo-tools

python subsamp_mip.py $SLURM_ARRAY_TASK_ID


exit 0
######################################################
######################################################
######################################################




#sbatch --array=0-30 ./run_subsamp_mip.sh
sbatch --array=0-`grep padlock samplemeta.csv | wc -l` ./run_subsamp_mip.sh
