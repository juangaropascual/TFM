import mdtraj as md
import numpy as np
import sys
import csv
import glob


#Recognise each argument
for i in range(0, len(sys.argv)):
    if sys. argv[i] == "--chainA":
        chainA = sys.argv[i + 1]
    elif sys.argv[i] == "--chainB":
        chainB = sys.argv[i + 1]
    elif sys.argv[i] == "--out":
        out_pdb = sys.argv[i + 1]
    elif sys.argv[i] == "--smiles":
        SMILES = sys.argv[i + 1]
    elif sys.argv[i] == "--replicas":
        replicas = int(sys.argv[i + 1])
    if len(sys.argv) < 9:
        print('Error: python auto_diffdock.py --chainA chainA.pdb --chainB chainB.pdb --out out_file.pdb --smiles SMILES --replicas n_replicas')
        exit()


################################################
#                                              #
#               MERGE PDB FILES                #
#                                              #
################################################

#Asegurarse que las estructuras estan separadas
def distance(xyz1, xyz2):

    '''
    This function ensures that the structures are separated by at least 10 angstroms.

    Inputs:
    xyz1 - The coordinates of the first structure.
    xyz2 - The coordinates of the second structure.

    Output:
    xyz2 - The coordinates of the second structure, shifted by 10 angstroms if necessary. 
      
    '''

    if np.linalg.norm(xyz1.mean(axis = 1) - xyz2.mean(axis = 1)) < 10:
        mean1 = xyz1.mean(axis = 1)
        mean2 = xyz2.mean(axis = 1)
        vector = mean2 - mean1
        distance = np.linalg.norm(vector)
        unit_vector = vector / distance
        separation_vector = unit_vector * 10
        xyz2 += separation_vector


#Load the structures
p1 = md.load(chainA) 
p2 = md.load(chainB)

top1, xyz1 = p1.top, p1.xyz
top2, xyz2 = p2.top, p2.xyz
top1_2 = top1.join(top2)

distance(xyz1, xyz2)

#Merge the structures and save the new pdb file
xyz1_2 = np.concatenate((xyz1, xyz2), axis=1)
p1_2 = md.Trajectory(xyz=xyz1_2, topology=top1_2) 
p1_2.save(out_pdb)
print('New pdb file created: ' + out_pdb)


################################################
#                                              #
#               CREATE .csv FILE               #
#                                              #
################################################


def create_csv(protein_path, ligand_description):

    '''
    Creates a .csv file with the protein path and SMILES.

    Imputs:
    protein_path - A string representing the path to the protein. (e.g. '1a07.pdb') si este archivo esta en pdb_files
    ligand_description - A string representing the SMILES of the ligand. (e.g. 'CNCC1=CC=C(C=C1)C2=C3CCNC(=O)C4=C3C(=CC(=C4)F)N2')

    Output:
    protein_ligand.csv - A .csv file with the columns: complex_name, protein_path, ligand_description, protein_sequence.   
    
    '''

    with open('protein_ligand.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
        for i in range(replicas):
            writer.writerow(["", protein_path, ligand_description, ""])

file = out_pdb # es lo mismo que se pone en --out del merge.py pero con .pdb
#SMILES = 'CNCC1=CC=C(C=C1)C2=C3CCNC(=O)C4=C3C(=CC(=C4)F)N2'
create_csv(file, SMILES)
print('New .csv file created: protein_ligand.csv')

