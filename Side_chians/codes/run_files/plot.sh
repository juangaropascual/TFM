#!/bin/bash
source input.dat

python codes/plots.py --run_name "$run_name" --AF "$AF" --DF "$DF" --OF "$OF"