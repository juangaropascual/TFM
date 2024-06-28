#!/bin/bash

source input.dat # Source the input variables

echo ''
pred_list=$(ls output | grep pred_*)

if [ -z "$pred_list" ]; then
        echo "There is no file starting with 'pred_' on the output directory."
        echo " --> Use 'make rank' to generate the necessary files."
        echo " --> Manually provide the files in the output directory."
fi
# plot the rest of the plots
INDEX=0
for pred in $pred_list; do
    if [ $INDEX == 0 ]; then
        python codes/plots.py --data_a "$data_A" --data_b "$data_B" --threshold "$th" --predictions output/"$pred" --failed_file "$failed" --plot_th $plot_th 
    else
        python codes/plots.py --data_a "$data_A" --data_b "$data_B" --threshold "$th" --predictions output/"$pred" --failed_file "$failed" --plot_th $plot_th
    fi
    let INDEX=${INDEX}+1
done