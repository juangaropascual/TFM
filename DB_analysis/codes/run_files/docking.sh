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
export PYTHONPATH=$DD_path #path to where DiffDock is installed
dir=$(pwd)
dd_dir=$DD_path
cd $dd_dir

if [ $docking_program == 'DD' ]; then

    echo "Starting Diffdock..."
    date

    python -W ignore inference.py \
    --protein_ligand_csv $dir/output/protein_ligand_$run_name.csv \
    --out_dir $dir/output/results_dd/$run_name \
    --inference_steps 20 \
    --samples_per_complex $samples \
    --batch_size 10 \
    --actual_steps 18 \
    --no_final_step_noise \
    --model_dir workdir/paper_score_model \
    --confidence_model_dir workdir/paper_confidence_model
fi

if [ $docking_program == 'DD_L' ]; then

    echo "Starting Diffdock-L..."
    date

    python -W ignore inference.py \
    --config default_inference_args.yaml \
    --protein_ligand_csv $dir/output/protein_ligand_$run_name.csv \
    --out_dir $dir/output/results_dd/$run_name 
fi

if [ $docking_program != 'DD' ] && [ $docking_program != 'DD_L' ]; then
    echo "Error: (input.dat) The variable docking_program must be 'DD' or 'DD_L'."
fi

cd $dir

if [ $auto_analysis == true ]; then
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
fi