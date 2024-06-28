#!/bin/bash

source input.dat # Source the input variables

# Clean database

if [ $new == false ]; then
    python codes/DB_cleaning.py --affinity "$affinity" --sequence "$sequence" --CID "$CID" --SMILES "$SMILES" --threshold "$threshold" --DB_file "$database" --match_file "$match_file"
else
    python codes/DB_cleaning.py --display_col_names True --affinity "$affinity" --sequence "$sequence" --CID "$CID" --SMILES "$SMILES" --threshold "$threshold" --DB_file "$database" --match_file "$match_file"
fi