#!/bin/bash

#SBATCH -A im3
#SBATCH -N 4
#SBATCH --ntasks-per-node 1
#SBATCH -p slurm
#SBATCH -t 05:00:00
#SBATCH --job-name farmbulk

module purge
module load python/miniconda4.12
source /share/apps/python/miniconda4.12/etc/profile.d/conda.sh
conda activate farmwell

srun --wait=0 driver.sh

conda deactivate
