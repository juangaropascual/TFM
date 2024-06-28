#!/bin/bash
source codes/run_files/env.sh # Activate env if necessary
source input.dat #source the input variables   

if [ -z $SMILES ]; then 
    echo "Use input.dat to introduce the SMILES of the ligand to start."
fi

job_names=$(ls protein_structures)

# if [[ ! "${job_names[@]}" =~ "$run_name" ]]; then
#     # random_int=$((RANDOM % 9000 + 1000))
#     # run_name="job_$random_int"
#     # sed -i "s/^run_name=.*/run_name=\"$run_name\"/" input.dat
#     # echo "Assigning random name to job since none was provided in input.dat."
#     echo "run_name must correspond to one of the protein_structures directories: "
#     echo "$job_names"
# else
#     echo ""
#     echo ""

#     echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
#     echo "                $run_name"
#     echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
#     echo ""
#     alias peligand=make
#     make $@
# fi

echo ""
echo ""

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "                $run_name"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo ""
alias peligand=make
make $@