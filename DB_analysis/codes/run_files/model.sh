#!/bin/bash

source input.dat # Source the input variables

echo ''
results_dd=$(ls output/results_dd)

if [ ! -z "$results_dd" ]; then
        echo "There are no DiffDock results in results_dd."
        echo " --> Use 'diffduck run_dd' to generate the necessary files."
        echo " --> Manually provide the files in the output directory."
else
    previous_labels=$(ls -p model/training_data | grep -v /)
    if [ $new_labels == true ] || [ -z "$previous_labels" ]; then
        echo "New labels will be generated."
        python codes/mix.py --data_a "$data_A" --data_b "$data_B" --th_amb "$th_amb" --labels "$labels"
    fi

    for file in $(ls -p model/training_data | grep -v /); do
        python codes/prepare_features.py --training_data model/training_data/"$file" --n_samples "$n_samples" --labels "$labels"
    done
    
    echo 'Starting forest training... '
    for dir in $(ls -p model/training_data | grep /); do
        python codes/generate_forest_model.py --data model/training_data/"$dir" --n_samples "$n_samples"  --model_name "$model_name" --labels "$labels"
        # python codes/svm_model.py --data model/training_data/"$dir" --n_samples "$n_samples"  --model_name "$model_name"
    done
fi

