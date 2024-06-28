#!/bin/bash

source input.dat # Source the input variables

result_list=$(ls DD)

echo ''

if [ -z "$result_list" ]; then
    echo "/DD is empty."
    echo " --> Use output/Af_input.csv to run the predicted structures for the proteins"
    echo " --> Use the predicted proteins to generate the DD predictions."
    echo " --> If you already have a results file (i.e results_DD_AF.csv) place it in the /DD directory."
    echo " --> When there is at least one file in the /DD directory the code will generate the ranking and the plots."

fi

# Rank the predictions and plot ki ratio vs rel
for result in $result_list; do
    python codes/rank_formula.py --results_file DD/"$result" --output_file pred_"$result" --data_a "$data_A" --data_b "$data_B"
done