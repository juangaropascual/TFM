#!/bin/bash
source input.dat # Source the input variables

sbatch -J $run_name codes/run_files/analysis.sh