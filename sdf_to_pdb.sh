#!/bin/bash
#path="/home/ramon/juan/diffdock_test/results/prova4.0/ligands"
path=$1

for file in $path/*sdf
do
    base=$(basename "$file" .sdf)
    obabel "$file" -O "$path/${base}.pdb" 
done
echo "Conversion finished"
rm $path/*.sdf
echo "SDF files removed"


