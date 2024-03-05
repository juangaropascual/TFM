#!/bin/bash
# ./run_analysis.sh /folder_of_ligands protein.pdb

folder=$1
pdb_file=$2
path_ligand=$PWD/$1
path=$PWD
cd /home/ramon/juan/important_files
./sdf_to_pdb.sh $path_ligand
python /home/ramon/juan/important_files/ranking_afinity.py $path/$pdb_file $path_ligand


