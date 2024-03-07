# Diffdock: run and analysis of binding of one ligand in a complex with two proteins 

With this implementation of Diffdock the relative affinity of two molecules to the same legand can be analysed. The arguments needed to run the analysis `all_run.sh` are `--chainA`, a PDB file that will be analised as the chain A, `--chainB`, a PDB file that will be analysed as chain B. `--out`, a PDB file that will contain the chain A and B merged in one file, the name of the file will be the name of the folder with the results of analysis. `--smiles`, the SMILES of the ligand and `--replicas` the number of times you whany to run the Diffdock and the analysis. 

## Usage

```
./all_run.sh --chainA chainA.pdb --chainB chainB.pdb --out output.pdb --smiles 'SMILES' --replicas n
```

All the PDB files must be "clean", without water molecules, and they must be only one chain per PDB file for the analysis to work.
