#!/bin/bash
source input.dat
echo "$folding"
python codes/rank.py --run_name "$run_name" --AF "$AF" --DF "$DF" --OF "$OF" --n_samples "$n_samples"