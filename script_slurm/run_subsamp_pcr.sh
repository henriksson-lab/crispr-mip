#!/usr/bin/env bash
#SBATCH -A naiss2024-22-1647
#SBATCH -J pcr
#SBATCH -n 64
#SBATCH -t 0-24:00:00
#SBATCH -p shared
#SBATCH --mem=100240

######################################################
#module load bioinfo-tools

python subsamp_pcr.py $SLURM_ARRAY_TASK_ID


exit 0
######################################################
######################################################
######################################################




#lib pcr
#sbatch --array=0-5 ./run_subsamp_pcr.sh

sbatch --array=0-`grep PCR samplemeta.csv | wc -l` ./run_subsamp_pcr.sh
