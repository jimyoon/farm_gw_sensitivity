
# creating the environment:
- `module load python/miniconda4.12`
- `source /share/apps/python/miniconda4.12/etc/profile.d/conda.sh`
- `conda create -n farmwell python=3.11`
- `conda activate farmwell`
- `conda install -c conda-forge numpy pandas scipy SALib pyomo matplotlib openpyxl ipopt`
- `conda deactivate`

# running the experiment
- Update the total number of nodes to use in `q.sl`
- Update the expected walltime in `q.sl`
- Update the partition in `q.sl`
- Update the total number of farms to run in `driver.sh`
- `sbatch q.sl`
- results and errors are written to `./output` directory per farm
