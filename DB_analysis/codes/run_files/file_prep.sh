#!/bin/bash

# source codes/run_files/env.sh # Activate env if necessary
source input.dat # Source the input variables

python codes/input_merge.py --chain_A $data_A --chain_B $data_B --name $run_name
echo "--> Merge input file created"

mkdir -p output/prot_structures/merged_structures/$run_name
echo "--> Merged structures folder created ($run_name)"
echo "--> Merging PDB files..."

python codes/file_preparation_multi.py --merge_input "output/merge_input.csv" --replicas $n_replica --name $run_name
echo "--> protein_ligand_$run_name.csv created (Input for DiffDock)"
echo "--> File preparation done"