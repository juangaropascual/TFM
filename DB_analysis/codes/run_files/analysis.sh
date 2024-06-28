#!/bin/bash

source input.dat # Source the input variables

analysis_counter=0
result_dir=output/results_dd/$run_name/
structure_dir=output/prot_structures/merged_structures/
for folder in $result_dir*; do
    if [ -d "$folder" ]; then
        echo "Analysing "$folder""
        file_name=$(basename "$folder")
        codes/sdf_to_pdb.sh "$folder"
        python codes/ranking_afinity_new.py --protein "${structure_dir}${run_name}/${file_name}.pdb" --folder "$folder" --counter $analysis_counter
        analysis_counter=$((analysis_counter+1))
    fi

echo "-->Done "$folder""
echo "-->Results are in $"$file_name""
echo "-------------------------------------------------------------------"
done

cd output/results_dd/${run_name}
awk 'FNR==1 && NR!=1{next;}{print}' $(find . -name 'result.csv') > ../../../DD/result_${run_name}.csv
echo " --> Analysis finished"
echo " --> Results are in DD/result_{$run_name}.csv"
