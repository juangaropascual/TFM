#!/bin/bash
source input.dat # Source the input variables

mkdir -p output/results_dd/$run_name
sbatch -J $run_name codes/run_files/docking.sh