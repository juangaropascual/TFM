#!/bin/bash

chainA=$1
chainA_file=$2
chainB=$3
chainB_file=$4
out=$5
out_file=$6
smiles=$7
smiles_string=$8
replicas=$9
n_replicas=${10}



python /home/ramon/juan/important_files/file_preparation.py $chainA $chainA_file $chainB $chainB_file $out $out_file $smiles $smiles_string $replicas $n_replicas
echo "Starting Diffdock..."

exec &> analysis.log

python -W ignore /home/ramon/progs/DiffDock/inference.py \
--protein_ligand_csv protein_ligand.csv \
--out_dir results/ \
--inference_steps 20 \
--samples_per_complex 40 \
--batch_size 10 \
--actual_steps 18 \
--no_final_step_noise \
--model_dir /home/ramon/progs/DiffDock/workdir/paper_score_model \
--confidence_model_dir /home/ramon/progs/DiffDock/workdir/paper_confidence_model 

analysis_counter=0
result_dir=results/
for folder in $result_dir*; do
    if [ -d "$folder" ]; then
        echo "$folder"
        /home/ramon/juan/important_files/sdf_to_pdb.sh "$folder"
        python /home/ramon/juan/important_files/ranking_afinity_def.py $out_file "$folder" $analysis_counter
        analysis_counter=$((analysis_counter+1))
    fi
done

exec &> /dev/tty
echo "Done"
mkdir "${out_file%.*}"
mv results/complex_* "${out_file%.*}"
mv *csv *log $out_file "${out_file%.*}"
mv .log 6xzs-1lzv/
echo "Results are in ${out_file%.*}"