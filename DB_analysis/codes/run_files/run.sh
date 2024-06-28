#!/bin/bash

# source codes/run_files/env.sh # Activate env if necessary
source input.dat # Source the input variables

# Clean database

if [ $new == false ]; then
    python codes/DB_cleaning.py --affinity "$affinity" --sequence "$sequence" --CID "$CID" --SMILES "$SMILES" --threshold "$threshold" --DB_file "$database" --match_file "$match_file"
else
    python codes/DB_cleaning.py --display_col_names True --affinity "$affinity" --sequence "$sequence" --CID "$CID" --SMILES "$SMILES" --threshold "$threshold" --DB_file "$database" --match_file "$match_file"
fi


if [ $new == false ]; then
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

    echo ''

    pred_list=$(ls output | grep pred_*)

    # plot the rest of the plots
    INDEX=0
    for pred in $pred_list; do
        if [ $INDEX == 0 ]; then
            python codes/plots.py --data_a "$data_A" --data_b "$data_B" --threshold "$th" --predictions output/"$pred" --counts True --failed_file "$failed" --plot_th $plot_th 
        else
            python codes/plots.py --data_a "$data_A" --data_b "$data_B" --threshold "$th" --predictions output/"$pred" --failed_file "$failed" --plot_th $plot_th
        fi
        let INDEX=${INDEX}+1
    done

    # Mix data A wit data B in order to prevent errors while training
    python codes/mix.py --data_a "$data_A" --data_b "$data_B" --th_amb "$th_amb"
fi

