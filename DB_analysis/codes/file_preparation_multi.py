import mdtraj as md
import numpy as np
from argparse import ArgumentParser
import pandas as pd
import os


parser = ArgumentParser()
parser.add_argument(
    "--merge_input",
    type=str,
    help="CSV file with the protein-ligand information. The columns are: complex_name, chainA, chainB, SMILES. (output/merge_input.csv)",
    required=True,
    default='output/merge_input.csv'
)
parser.add_argument(
    "--replicas",
    type=str,
    help="Number of replicas to generate",
    required=False,
    default= 1
)
parser.add_argument(
    "--name",
    type=str,
    help="Name of the run",
    required=True,
)
args = parser.parse_args()

################################################
#                                              #
#               MERGE PDB FILES                # 
#                                              #
################################################
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


################################################
#                                              #
#               CREATE .csv FILE               # 
#                                              #
################################################
def create_csv(data, protein_complex):

    '''
    Creates a .csv file with the protein path and SMILES, that will be used as input for Diffdock.

    Inputs:
    data - A .csv file with the columns: complex_name, chainA, chainB, SMILES.

    Output:
    protein_ligand.csv - A .csv file with the columns: complex_name, protein_path, ligand_description, protein_sequence.   
    
    '''

    protein_ligand = {'complex_name':[], 'protein_path': [], 'ligand_description': [], 'protein_sequence': []}
    protein_ligand['complex_name'] = protein_complex
    protein_ligand['protein_path'] = current_path + '/' + path + protein_complex + '.pdb'
    protein_ligand['ligand_description'] = data['SMILES']
    protein_ligand['protein_sequence'] = np.nan

    df = pd.DataFrame(protein_ligand)
    df.to_csv('output/protein_ligand_' + run_name + '.csv', index=False)

# Main code
if __name__ == '__main__':
    merge_input = args.merge_input
    replicas = args.replicas
    run_name = args.name
    path = 'output/prot_structures/merged_structures/' + run_name + '/'
    current_path = os.getcwd()
    data = pd.read_csv(merge_input)
    protein_complex = data['complex_name']
    data_array = data.values


    for complex_name, chainA, chainB, SMILES in data_array:

        p1 = md.load(chainA)
        p2 = md.load(chainB)

        top1, xyz1 = p1.top, p1.xyz
        top2, xyz2 = p2.top, p2.xyz
        top1_2 = top1.join(top2)

        distance(xyz1, xyz2)

        xyz1_2 = np.concatenate((xyz1, xyz2), axis=1)
        p1_2 = md.Trajectory(xyz=xyz1_2, topology=top1_2) 
        p1_2.save(path + complex_name + '.pdb')
        # print('New pdb file created: ' + complex_name + '.pdb')
    create_csv(data, protein_complex)
   

