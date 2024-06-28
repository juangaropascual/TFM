#!/bin/bash
#SBATCH -J $run_name
#SBATCH -e error_%j.err
#SBATCH -o output_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=long
#SBATCH --gres=gpu:1
##SBATCH --constraint=gpu8

source input.dat # Source the input variables

dir=$(pwd)

input="$dir/inputs/folding_${run_name}_input.csv"
cd protein_structures/
mkdir -p ${run_name}
cd ${run_name}

#########################################
#              ALPHAFOLD                #
#########################################

if [ $AF == 1 ]; then

    mkdir -p AF
    mkdir -p AF/AF_res

    colabfold_batch $input "$dir"/protein_structures/"$run_name"/AF

    cd AF && $(cp `ls | grep rank_001 | grep .pdb | xargs` AF_res)
    cd AF_res
    for file in $(ls); do
        if [[ -f "$file" ]]; then 
            filename=$(basename "$file")
            new_filename=$(echo "$filename" | cut -d '_' -f 1)
            mv "$file" "$new_filename.pdb" 
        fi
    done
    cd ..
    rm -f *. 
    rm -r $run_name*
    mv AF_res/* .
    rm -r AF_res
    rm -f cite.bibtex log.txt config.json
    cd ..
fi

#########################################
#               DEEPFOLD                #
#########################################

if [ $DF == 1 ]; then
    mkdir -p DF
    mkdir -p DF/DF_res
    colabfold_batch --model-type deepfold_v1 $input "$dir"/protein_structures/"$run_name"/DF

    cd DF && $(cp `ls | grep rank_001 | grep .pdb | xargs` DF_res)
    cd DF_res
    for file in $(ls); do
        if [[ -f "$file" ]]; then 
            filename=$(basename "$file")
            new_filename=$(echo "$filename" | cut -d '_' -f 1)
            mv "$file" "$new_filename.pdb" 
        fi
    done
    cd ..
    rm -f *. 
    rm -r $run_name*
    mv DF_res/* .
    rm -r DF_res
    rm -f cite.bibtex log.txt config.json
    cd ..
fi

#########################################
#              OMEGAFOLD                #
#########################################

if [ $OF == 1 ]; then
    mkdir -p OF

    # Check if the last line is empty
    if [[ -z $(tail -n 1 "$input") ]]
    then
        sed -e :a -e '/^\n*$/{$d;N;};/\n$/ba' "$input" > temp.txt
        mv temp.txt "$input"
    else
        echo "" >> "$input"
    fi

    while IFS=, read -r sequence id Affinity
    do
        if [[ $id != "id" ]]
        then
            echo ">${id}"
            echo "${sequence}"
        fi
    done < $input > input.fasta

    omegafold input.fasta "$dir"/protein_structures/"$run_name"/OF
    rm input.fasta
fi
