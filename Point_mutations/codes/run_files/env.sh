#!/bin/bash

directory=$(pwd)

# Activate the conda environment for DD if necessary

if [[ $(conda env list | grep "*" | awk '{print$1}') != 'diffdock' ]]; then
    echo "Activating DiffDock conda environment"
	conda activate diffdock
fi
echo ''
echo '----------------------------------------------------------------'
# create necessary directories

# Check if /output directory exists
if [ -d "protein_structures" ]; then
    # If it exists, erase its contents
    echo -e  "The \e]8;;file:$directory/protein_structures\a protein_structures \e]8;;\a directory contains:"
    echo "$(ls -hs protein_structures)"
    # rm -r output/*
    # mkdir output/figures
else
    # If it doesn't exist, create it
    mkdir protein_structures
fi

echo '----------------------------------------------------------------'
echo ''
echo '----------------------------------------------------------------'

if [ -d "inputs" ]; then
    echo -e  "The \e]8;;file:$directory/inputs\a inputs \e]8;;\a directory contains:"
    echo "$(ls -hs inputs)"
else
    # If it doesn't exist, create it
    mkdir inputs
fi
echo '----------------------------------------------------------------'
echo ''
echo '----------------------------------------------------------------'

if [ -d "results" ]; then
    echo -e  "The \e]8;;file:$directory/results\a results \e]8;;\a directory contains:"
    echo "$(ls -hs results)"
else
    # If it doesn't exist, create it
    mkdir results
    mkdir results/figures
fi
echo '----------------------------------------------------------------'

