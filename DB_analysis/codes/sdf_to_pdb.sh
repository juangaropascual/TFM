#!/bin/bash
path=$1

for file in $path/*sdf
do
    base=$(basename "$file" .sdf)
    obabel "$file" -O "$path/${base}.pdb" 
done
echo "Conversion finished"
rm $path/*.sdf
rm $path/rank1.pdb
echo "SDF files removed"