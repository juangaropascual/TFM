#!/bin/bash

directory=$(pwd)

# Activate the conda environment for DD if necessary

if [[ $(conda env list | grep "*" | awk '{print$1}') != 'diffdock' ]]; then
    echo "Activating conda environment"
	conda activate diffdock
fi
echo ''
echo '----------------------------------------------------------------'
# create necessary directories

# Check if /output directory exists
if [ -d "output" ]; then
    # If it exists, erase its contents
    echo -e  "The \e]8;;file:$directory/output\a output \e]8;;\a directory already exists:"
    echo "$(ls -hs output)" 
    # rm -r output/*
    # mkdir output/figures
else
    # If it doesn't exist, create it
    mkdir output
    mkdir output/figures
    mkdir output/results_dd
    mkdir output/prot_structures
    mkdir output/prot_structures/merged_structures
fi

echo '----------------------------------------------------------------'
echo ''
echo '----------------------------------------------------------------'
# Check if /DD directory exists
if [ ! -d "DD" ]; then
    # If it doesn't exist, create it
    mkdir DD
else
    echo -e  "The \e]8;;file:$directory/DD\a DD \e]8;;\a directory already exists:"
    echo "$(ls -hs DD)"
fi
echo '----------------------------------------------------------------'
echo ''
echo '----------------------------------------------------------------'
# Check if /model directory exists
if [ ! -d "model" ]; then
    # If it doesn't exist, create it
    mkdir model
    mkdir model/training_data
else
    echo -e  "The \e]8;;file:$directory/model\a model \e]8;;\a directory already exists:"
    echo "$(ls -hs model)"
fi
echo '----------------------------------------------------------------'
echo ''

# Iterate over each zip file found and unzip it
for zip_file in $(find "$directory" -type f -name "*.zip")
do
    # Unzip the file
    unzip "$zip_file" -d "${zip_file%.*}"
    echo "Unzipped: $zip_file"
done
