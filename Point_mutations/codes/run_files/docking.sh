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
mkdir -p results/results_dd/$run_name

fold=$1

if [ $docking_program == 'DD' ]; then

    echo "Starting Diffdock..."
    date

    python -W ignore inference.py \
    --protein_ligand_csv $dir/inputs/docking_${run_name}_"$fold"_input.csv \
    --out_dir $dir/results/results_dd/"$run_name"_"$fold" \
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

    python -W ignore -m inference \
    --config default_inference_args.yaml \
    --protein_ligand_csv $dir/inputs/docking_${run_name}_"$fold"_input.csv \
    --out_dir $dir/results/results_dd/"$run_name"_"$fold" --samples_per_complex $samples
fi

if [ $docking_program != 'DD' ] && [ $docking_program != 'DD_L' ]; then
    echo "Error: (input.dat) The variable docking_program must be 'DD' or 'DD_L'."
fi

cd $dir

if [ $auto_analysis == true ]; then
    echo "Analysing $run_name"
    analysis_counter=0
    result_dir=results/results_dd/"$run_name"_"$fold"/
    structure_dir=protein_structures/$run_name/merged_structures/$fold 
    for folder in $result_dir*; do
        if [ -d "$folder" ]; then
            echo "Analysing "$folder""
            file_name=$(basename "$folder")
            codes/sdf_to_pdb.sh "$folder"
            python codes/ranking_afinity_new.py --protein "${structure_dir}/${file_name}.pdb" --folder "$folder" --counter $analysis_counter
            analysis_counter=$((analysis_counter+1))
        fi

    echo "--> Done $folder"
    echo "--> Results are in $file_name"
    echo "-------------------------------------------------------------------"
    done

    cd results/results_dd/${run_name}_$fold
    awk 'FNR==1 && NR!=1{next;}{print}' $(find . -name 'result.csv') > ../../../results/docking_${run_name}_$fold.csv
    echo " --> Analysis finished"
    echo " --> Results are in results/result_${run_name}_$fold.csv"
fi