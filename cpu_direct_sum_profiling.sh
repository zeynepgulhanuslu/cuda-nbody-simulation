#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 output_directory"
    exit 1
fi

OUTDIR=$1
mkdir -p "$OUTDIR"

N_LIST=(200 300 500 1000 2000 100000 500000)
STEPS_LIST=(1000 2000 3000 5000)
RESULT_FILE="${OUTDIR}/all_result.txt"
: > "$RESULT_FILE"

for N in "${N_LIST[@]}"; do
    for STEPS in "${STEPS_LIST[@]}"; do
        OUTFILE="${OUTDIR}/N${N}_steps${STEPS}.hdf5"
        PROFILE_OUT="${OUTDIR}/cpu_profile_N${N}_steps${STEPS}.txt"
        echo "Running: N=$N, steps=$STEPS"

        run_output=$(./nbody_cpu "$OUTFILE" "$N" "$STEPS")
        total_time=$(echo "$run_output" | grep "Process completed! Total time:" | awk '{print $(NF-1)}')

        echo "N=$N, Step=$STEPS, total_time=$total_time" | tee -a "$RESULT_FILE"
        gprof ./nbody_cpu gmon.out > "$PROFILE_OUT"
        rm -f gmon.out
    done
done
