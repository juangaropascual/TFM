#!/bin/bash
source input.dat # Source the input variables

clust=$(which sbatch) # check if cluster

if [ -z "$clust" ]; then

    if [ $AF == 1 ]; then
        ./codes/run_files/structures.sh
    fi

else

    if [ $AF == 1 ]; then
        sbatch -J fold_"$run_name" codes/run_files/structures.sh
    fi

fi