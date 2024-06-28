#!/bin/bash
source input.dat

python codes/mutate_DS.py --run_name "$run_name" --ds "$ds" --n_mutations "$n_mutations"