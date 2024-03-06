#!/bin/bash

#SBATCH -A im3
#SBATCH -N 2
#SBATCH --ntasks-per-node 1
#SBATCH -p short
#SBATCH -t 00:30:00
#SBATCH --job-name farmwell

module purge
module load python/miniconda4.12
source /share/apps/python/miniconda4.12/etc/profile.d/conda.sh
conda activate farmwell

srun --wait=0 driver.sh

conda deactivate
