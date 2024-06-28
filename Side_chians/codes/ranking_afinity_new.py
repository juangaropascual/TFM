import mdtraj as md
import os
from itertools import product
import numpy as np
import re
import pandas as pd
from argparse import ArgumentParser
import glob

parser = ArgumentParser()
parser.add_argument(
    "--protein",
    type=str,
    help="PDB file of the protein complex (protein_complex.pdb)",
    required=True,
)
parser.add_argument(
    "--folder",
    type=str,
    help="Directory of the Diffdock results where all the ligand predictions are for each protein complex (protein_complex)",
    required=True,
)
parser.add_argument(
    "--counter",
    type=str,
    help="Number of folders analysed. Each directory analysed increase the counter",
    required=True,
)
args = parser.parse_args()

def load_ligand(path, file):

    """
    This function loads a structure file, typically representing the ligand.
    
    Inputs:
    path - A string that represents the directory where the file is located.
    file - A string that represents the name of the file.
    
    The function concatenates `path` and `file` with a '/' in between to create the full path to the file. 
    This full path is then passed to the `md.load` function from the `mdtraj` library.
    
    Output:
    The loaded ligand.
    """

    ligand = md.load(path + '/' + file)
    return ligand

def extract_number(filename):

    """
    This function extracts the first sequence of digits from a given filename. To sort the ligands in numeric order.
    
    Input:
    filename - A string that represents the name of the file.
    
    The function uses the `re.search` function from Python's `re` module to find one or more digits in the `filename`.
    If a match is found, the `match.group()` function is used to return the actual matched text, which is then converted to an integer.
    
    Output:
    Order of ligands, or 0 if no sequence of digits is found.
    """

    match = re.search(r'\d+', filename)
    return int(match.group()) if match else 0

def extract_confidence(filename):

    """
    This function extracts the confidence value from the Diffdock predictions.
    
    Input:
    filename - A string that represents the name of the file.
    
    The function uses the `re.search` function from Python's `re` module to find a sequence that starts with 'confidence', 
    followed by one or more digits or a period, which can optionally start with a minus sign, in the `filename`.
    If a match is found, the `match.group(1)` function is used to return the actual matched text, which is then returned as a string.
    
    Output:
    The confidence value as a string, or 0 if no confidence value is found.
    """

    match = re.search(r'confidence(-?[\d]+.[\d][\d])', filename)
    return str(match.group(1)) if match else 0

def calculate_distances(prot, ligand):
    """
    This function calculates the minimum distances between a protein and a ligand.

    Parameters:
    prot (md.Trajectory): The protein as a molecular dynamics trajectory.
    ligand (md.Trajectory): The ligand as a molecular dynamics trajectory.

    The function performs the following tasks:
    1. It extracts the topology and coordinates from the protein and ligand trajectories.
    2. It joins the topologies and concatenates the coordinates to create a new trajectory that includes both the protein and ligand.
    3. It selects the alpha carbons (CA) of chain A and chain B of the protein and all atoms of the ligand (residue name UNL).
    4. It creates pairs of alpha carbons and ligand atoms for chain A and chain B separately.
    5. It computes the distances between the pairs of atoms for chain A and chain B separately and finds the minimum distance in each case.

    Returns:
    dist_A (float): The minimum distance between the alpha carbons of chain A and the ligand.
    dist_B (float): The minimum distance between the alpha carbons of chain B and the ligand.
    """

    top_prot, xyz_prot = prot.top, prot.xyz
    top_lig, xyz_lig = ligand.top, ligand.xyz
    top_all = top_prot.join(top_lig)
    xyz_all = np.concatenate((xyz_prot, xyz_lig), axis=1)

    prot_lig = md.Trajectory(xyz=xyz_all, topology=top_all)

    chainA = prot_lig.top.select('chainid 0 name CA')
    chainB = prot_lig.top.select('chainid 1 name CA')
    ligand = prot_lig.top.select('resname UNK')
    if ligand.size == 0:
        ligand = prot_lig.top.select('resname UNL')    
    pairs_A = list(product(chainA, ligand))
    dist_A = md.compute_distances(prot_lig, pairs_A).min()

    pairs_B = list(product(chainB, ligand))
    dist_B = md.compute_distances(prot_lig, pairs_B).min()

    return dist_A, dist_B

def write_distances_to_csv(path, analysis_counter):
    """
    This function writes the data contained in the 'data_A' and 'data_B' variables to two separate CSV files.
    
    Parameters:
    path (str): The directory where the CSV files will be saved.
    analysis_counter (int): A counter used to differentiate the names of the CSV files.

    The function performs the following tasks:
    1. It creates two pandas DataFrames, 'df_A' and 'df_B', from the data in the 'data_A' and 'data_B' variables respectively. 
       'data_A' and 'data_B' are expected to be dictionaries or lists of lists (or similar data structures) 
       that contain the data to be written to the CSV files.
    2. It writes the 'df_A' DataFrame to a CSV file named 'distancias_A_{analysis_counter}.csv' and the 'df_B' DataFrame 
       to a CSV file named 'distancias_B_{analysis_counter}.csv' in the directory specified by the 'path' variable. 
       The 'analysis_counter' is used to differentiate the names of the CSV files. The mode 'a' is used to append data 
       if the file already exists. The index of the DataFrame is not written into the file.

    Note: 'data_A' and 'data_B' are not passed as arguments to the function, so they must be defined in the same 
    scope as this function for it to work correctly.
    """

    df_A = pd.DataFrame(data_A)
    df_B = pd.DataFrame(data_B)

    df_A.to_csv(f'{path}/distancias_A_{analysis_counter}.csv', mode='a', index=False)
    df_B.to_csv(f'{path}/distancias_B_{analysis_counter}.csv', mode='a', index=False)

def analyze_distances(distancias_A, distancias_B):

    """
    This function analyzes the distances between a protein and a ligand for two chains (A and B), and prints out the results.
    
    Inputs:
    distancias_A - A list of tuples, where each tuple contains the index, distance, and confidence for chain A.
    distancias_B - A list of tuples, where each tuple contains the index, distance, and confidence for chain B.
    
    The function iterates over the tuples in `distancias_A` and `distancias_B`. If the index (step) in a tuple is 0, 
    it prints that the rank 1 ligand is in the corresponding chain.
    
    The function then compares the lengths of `distancias_A` and `distancias_B`. If one list is longer than the other, 
    it prints that the corresponding chain has more ligands and calculates the percentage of ligands in that chain. 
    If the lengths are equal, it prints that both chains have the same number of ligands.
    """

    for tupla in distancias_A:
        step, distance, confidence = tupla
        if step == 0:
            print('Rank 1 ligand is in chain A')
    if len(distancias_A) > len(distancias_B):
        print('Chain A has more ligands')
        porcentaje = (len(distancias_A)/(len(distancias_A)+len(distancias_B))*100)
        print(porcentaje , '%')
    for tupla in distancias_B:
        step, distance, confidence = tupla
        if step == 0:
            print('Rank 1 ligand is in chain B')
    if len(distancias_B) > len(distancias_A):
        print('Chain B has more ligands')
        porcentaje = (len(distancias_B)/(len(distancias_A)+len(distancias_B))*100)
        print(porcentaje , '%')
    if len(distancias_B) == len(distancias_A):
        print('Both chains have the same number of ligands')

def write_results_to_csv(path,data_result):
    """
    This function writes the data contained in the 'data_result' variable to a CSV file.
    
    Parameters:
    path (str): The directory where the result.csv file will be saved.

    The function performs two main tasks:
    1. It creates a pandas DataFrame 'df_result' from the data in the 'data_result' variable. 
       'data_result' is expected to be a dictionary or a list of lists (or a similar data structure) 
       that contains the data to be written to the CSV file.
    2. It writes the 'df_result' DataFrame to a CSV file named 'result.csv' in the directory specified 
       by the 'path' variable.

    Note: 'data_result' is not passed as an argument to the function, so it must be defined in the same 
    scope as this function for it to work correctly.
    """

    df_result = pd.DataFrame(data_result)
    df_result.to_csv(f'{path}/result.csv')
    
if __name__ == '__main__':
    prot = md.load(args.protein)
    name_prot = args.protein
    work_index= name_prot.split('/')[-1].split('.')[0]
    print('Protein:', work_index)
    path = args.folder
    data_A = {'Step': [], 'Distance': [], 'Confidence': []}
    data_B = {'Step': [], 'Distance': [], 'Confidence': []}
    data_result = {'ID': [], 'Confidence': [], 'Chain(A=0)(B=1)': []}

    analysis_counter = args.counter
    files = [f for f in os.listdir(path) if f.endswith('.pdb')]    
    sorted_files = sorted(files, key=extract_number)
    
    for i, file in enumerate(sorted_files):
        ligand = load_ligand(path, file)

        confidence = extract_confidence(file)
        dist_A, dist_B = calculate_distances(prot, ligand)
        if dist_A < 0.4:
            data_A['Step'].append(i)
            data_A['Distance'].append(dist_A)
            data_A['Confidence'].append(confidence)
        if dist_B < 0.4:
            data_B['Step'].append(i)
            data_B['Distance'].append(dist_B)
            data_B['Confidence'].append(confidence)
        data_result['ID'].append(work_index)
        data_result['Confidence'].append(confidence)
        if dist_A < 0.4:
            data_result['Chain(A=0)(B=1)'].append(0)
        elif dist_B < 0.4:
            data_result['Chain(A=0)(B=1)'].append(1)
        else:
            data_result['Chain(A=0)(B=1)'].append(np.nan)
    write_distances_to_csv(path, analysis_counter)
    write_results_to_csv(path,data_result)
