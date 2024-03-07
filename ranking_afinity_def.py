import mdtraj as md
import os
from itertools import product
import numpy as np
import re
import sys
import csv

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

def calculate_distances(prot, ligand, i, confidence, distancias_A, distancias_B):

    """
    This function calculates the minimum distances between a protein and a ligand, and appends the index, distance, 
    and confidence to a list if the distance between ligand and protein is less than 0.4 Å.
    
    Inputs:
    prot - A protein object that contains topology and coordinates.
    ligand - A ligand object that contains topology and coordinates.
    i - An integer representing the index.
    confidence - A value representing the confidence.
    distancias_A - A list to store the index, distance, and confidence for chain A if the distance is less than 0.4.
    distancias_B - A list to store the index, distance, and confidence for chain B if the distance is less than 0.4.
    
    The function first extracts the topology and coordinates from the protein and ligand. It then joins the topologies 
    and concatenates the coordinates to create a new trajectory that includes both the protein and ligand.
    
    It selects the alpha carbons (CA) from chain A and chain B, and the ligand (resname UNL). It then calculates the 
    pairs of alpha carbons and ligand atoms for both chains.
    
    The function calculates the minimum distances between the pairs for both chains. If the distance for a chain is less 
    than 0.4, it appends a tuple of the index, distance, and confidence to the corresponding list.

    Output:
    Generates two csv file with the distances and confidence values.

    """

    top_prot, xyz_prot = prot.top, prot.xyz
    top_lig, xyz_lig = ligand.top, ligand.xyz
    top_all = top_prot.join(top_lig)
    xyz_all = np.concatenate((xyz_prot, xyz_lig), axis=1)

    prot_lig = md.Trajectory(xyz=xyz_all, topology=top_all)

    chainA = prot_lig.top.select('chainid 0 name CA')
    chainB = prot_lig.top.select('chainid 1 name CA')
    ligand = prot_lig.top.select('resname UNL')

    dist_A_file = f'distancias_A_{analysis_counter}.csv'
    pairs_A = list(product(chainA, ligand))
    dist_A = md.compute_distances(prot_lig, pairs_A).min()
    if dist_A < 0.4:
        distancias_A.append((i,dist_A,confidence))
        with open(dist_A_file, 'a', newline='') as f:  # Abre el archivo en modo de escritura
            writer = csv.writer(f)
            if os.stat(dist_A_file).st_size == 0:  # Si el archivo está vacío, escribe los encabezados
                writer.writerow(['Step', 'Distance', 'Confidence'])
            writer.writerow((i,dist_A,confidence))  # Escribe la tupla en el archivo

    dist_B_file = f'distancias_B_{analysis_counter}.csv'
    pairs_B = list(product(chainB, ligand))
    dist_B = md.compute_distances(prot_lig, pairs_B).min()
    if dist_B < 0.4:
        distancias_B.append((i,dist_B,confidence))
        with open(dist_B_file, 'a', newline='') as f:  # Abre el archivo en modo de escritura
            writer = csv.writer(f)
            if os.stat(dist_B_file).st_size == 0:  # Si el archivo está vacío, escribe los encabezados
                writer.writerow(['Step', 'Distance', 'Confidence'])
            writer.writerow((i,dist_B,confidence))  # Escribe la tupla en el archivo

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

    match = re.search(r'confidence(-?[\d.]+)', filename)
    return str(match.group(1)) if match else 0

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
    

# Main code
if __name__ == '__main__':
    prot = md.load(sys.argv[1])
    path = sys.argv[2]
    distancias_A = []
    distancias_B = []
    # i = 0
    analysis_counter = sys.argv[3]
    files = os.listdir(path)
    sorted_files = sorted(files, key=extract_number)

    for i,file in enumerate(sorted_files):
        ligand = load_ligand(path, file)
        confidence = extract_confidence(file)
        calculate_distances(prot, ligand, i, confidence, distancias_A, distancias_B)
        # i += 1

    analyze_distances(distancias_A, distancias_B)