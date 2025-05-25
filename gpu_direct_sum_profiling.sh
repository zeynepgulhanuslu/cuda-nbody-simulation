#!/bin/bash

# Usage: ./gpu_direct_sum_profiling.sh output_directory

if [ $# -ne 1 ]; then
    echo "Usage: $0 output_directory"
    exit 1
fi

OUTDIR=$1
mkdir -p "$OUTDIR"

N_LIST=(200 300 500 1000 2000 100000 500000)
STEPS_LIST=(1000 2000 3000 5000)

for N in "${N_LIST[@]}"; do
    for STEPS in "${STEPS_LIST[@]}"; do
        OUTFILE="${OUTDIR}/N${N}_steps${STEPS}_gpu.hdf5"
        PROFILE_OUT="${OUTDIR}/ncu_profile_N${N}_steps${STEPS}.ncu-rep"
        echo "Running (GPU): N=$N, steps=$STEPS"
        ncu --set full --target-processes all --export "$PROFILE_OUT" ./nbody_gpu "$OUTFILE" "$N" "$STEPS"
    done
done