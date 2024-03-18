#!/usr/bin/env bash
#SBATCH -A naiss2023-5-414
#SBATCH -J pcr
#SBATCH -n 16
#SBATCH -t 0-24:00:00
#SBATCH -p core

######################################################
#module load bioinfo-tools

python subsamp_pcr.py $SLURM_ARRAY_TASK_ID


exit 0
######################################################
######################################################
######################################################



sbatch --array=0-30 ./run_subsamp_pcr.sh




#lib pcr
sbatch --array=4-4 ./run_subsamp_pcr.sh

