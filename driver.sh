#!/usr/bin/env /bin/bash

# total number of farms to process
#total_farms=29679
total_farms=400

# total number of nodes used to run this job, defined in q.sl
total_nodes=$SLURM_JOB_NUM_NODES

start_farm=$((SLURM_NODEID*total_farms/total_nodes))
end_farm=$((SLURM_NODEID*total_farms/total_nodes+total_farms/total_nodes-1))

if (( end_farm >= total_farms )); then
  end_farm=$(total_farms-1)
fi

farms=( $(seq $start_farm $end_farm) )

echo "NODE $SLURM_NODEID: farms $start_farm - $end_farm"

export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2

/rcfs/projects/im3/gnuparallel/bin/parallel --jobs 32 "python run_experiment_HPC.py farm_ids_HPC_test.csv output/stephen_test_20240520 " ::: "${farms[@]}"

