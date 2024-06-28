#!/bin/bash

source input.dat # Source the input variables

# Clean database

python codes/create_input.py --SMILES $SMILES --run_name $run_name