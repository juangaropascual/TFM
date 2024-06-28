#!/bin/bash
source input.dat # Source the input variables

clust=$(which sbatch) # check if cluster

if [ -z "$clust" ]; then

    if [ $AF == 1 ]; then
        ./codes/run_files/docking.sh AF
    fi
    if [ $DF == 1 ]; then
        ./codes/run_files/docking.sh DF
    fi
    if [ $OF == 1 ]; then
        ./codes/run_files/docking.sh OF
    fi
else
    if [ $AF == 1 ]; then
        sbatch -J "$run_name"_AF codes/run_files/docking.sh AF
    fi
    if [ $DF == 1 ]; then
        sbatch -J "$run_name"_DF codes/run_files/docking.sh DF
    fi
    if [ $OF == 1 ]; then
        sbatch -J "$run_name"_OF codes/run_files/docking.sh OF
    fi
fi