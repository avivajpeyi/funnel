#!/bin/bash
#
#SBATCH --job-name=lartillot_1d
#SBATCH --output=logs/lartillot_1d_%A_%a.log
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=sstar
#SBATCH --array=0-5

module load {{module_loads}}
source {{python_env}}

ARRAY_ARGS=(0 1 2 3 4)

python run_lartillot.py -d 1 -s ${ARRAY_ARGS[${SLURM_ARRAY_TASK_ID}]}

