import pandas as pd
import numpy as np
from difflib import SequenceMatcher
from rdkit import Chem
import os

def read_db(file_path, columns, display_col_names=False):
    """
    Reads a tab-separated values file and returns a pandas dataframe.

    Args:
        file_path (str): The path to the file to be read.
        columns (List[str]): A list of the columns from the file that the user wants to extract.
        display_col_names (bool): Whether to display the column names. Default

    Returns:
        pandas.DataFrame: A pandas dataframe containing the selected data from the file.
    """
    if file_path.endswith(".tsv"):
        df = pd.read_csv(
            file_path, sep="\t", low_memory=False, on_bad_lines="skip"
        )  # Use '\t' as the delimiter
    elif file_path.endswith(".csv"):
        df = pd.read_csv(file_path, low_memory=False, on_bad_lines="skip")
    else:
        raise ValueError("File format not supported. Data base must be a tsv or csv")

    if display_col_names:
        for col in df.columns:
            if "Unnamed" not in col:
                print(col)
        return
    df = df[columns]
    return df

def seq_similarity(df):
    """
    This function compares the sequences of the different protein chains in the dataframe and prints the ones that are most similar.

    Args:
        df (pd.DataFrame): The dataframe containing the protein sequences and their affinities.

    Returns:
        None

    """

    global gl_seq
    global gl_affinity
    global gl_smiles

    red = '\033[31m'
    end = '\033[0m'

    un_lig = np.unique(list(df['SMILES']))
    lig_list = []
    input_list = []
    paths = []


    for lig in un_lig:
        sel = df[df['SMILES'] == lig]
        seq = list(sel[gl_seq])
        prot_names = list(sel['Target Name'])
        kd = list(sel[f"{gl_affinity} mean"])

        seq_list = []
        # aa_list = []
        kds = []
        for s1 in range(len(seq)):

            for s2 in range(s1 + 1, len(seq)):
                diff_count = 0
                # aa_mutated = 0
                seq1 = ''
                seq2 = ''
                # ind = 0

                if len(seq[s1]) == len(seq[s2]):
                
                    for a1,a2 in zip(seq[s1], seq[s2]):
                        # ind += 1
                        if a1!= a2:
                            diff_count += 1
                            seq1 += red+a1+end
                            seq2 += red+a2+end
                            # aa_mutated = ind
                        else: 
                            seq1 += a1
                            seq2 += a2
                    
                # diff_count = sum(a1 != a2 for a1, a2 in zip(seq[s1], seq[s2]))
                if diff_count == 1:
                    
                    if seq[s1] not in seq_list:
                        seq_list.append(seq[s1])
                        kds.append(kd[s1])
                        # aa_list.append(aa_mutated)

                    if seq[s2] not in seq_list:
                        seq_list.append(seq[s2])
                        kds.append(kd[s2])
                        # aa_list.append(aa_mutated)

                        print(prot_names[0].split(' ')[0], lig)
                        print(seq1)
                        print(seq2)
                        print('')
                    

        if len(seq_list) > 2:
            aux = pd.DataFrame()
            aux['sequence'] = seq_list
            
            # avoid different sequences in the same file.
            stop = False 
            for i in range(len(seq_list)):
                if len(seq_list[i]) != len(seq_list[0]):
                    stop = True 
            if stop:
                continue

            # Set an ID for every protein.
            id = []
            for i in range(len(seq_list)):
                if i == 0:
                    id.append(prot_names[0].split(' ')[0])
                else:
                    id.append(prot_names[0].split(' ')[0]+'-'+f'{i}')
            aux['id'] = id

            # Add the affinities
            aux[f"Affinity"] = kds 
            
            # Save the folding inputs in the inputs directory
            path = f"inputs/folding_{prot_names[0].split(' ')[0]}_input.csv"
            if path in paths:
                add = 1
                # If there are different smiles, add a number at the end.
                while path in paths: 
                    path = f"inputs/folding_{prot_names[0].split(' ')[0]}_input{add}.csv"
                    add +=1
                
            paths.append(path)
            aux.to_csv(path, index=False)

            input_list.append(path)
            lig_list.append(lig)

    # Create a file containing all the generated inputs and their smiles for the user to know which smiles corresponds to every folding input
    aux = pd.DataFrame()
    aux['Input'] = input_list
    aux['SMILES'] = lig_list
    aux.to_csv(f"inputs/search_smiles.csv", index=False)

def clean_df(
    df,
    min_atoms=8,
    max_atoms=60,
    seq_lim=500,
):
    """
    Cleans the dataframe by:
            --> Removing rows that have Nan values
            --> Removing rows that have Kd values including > symbol (i.e. >100000)
            --> Removing rows that contain ligands smaller than n_atoms
            --> Removing rows that contain Kd = 0
            --> Adding a column with the mean value of Kd for every lig-to-protein interaction
            --> Adding a column with the standard deviation of the mean Kd for every lig-to-protein interaction (if it is only one value then it is set to 0)
            --> Adding a column with the SMILES of every ligand
            --> Adding a column with the protein IDs

    Args:
        df (pandas.DataFrame): The dataframe to be cleaned.
        n_atoms (int, optional): The minimum number of atoms in a ligand for it to be included in the dataframe. Defaults to 8.
        seq_lim (int, optional): the maximum fasta length for a protein to be selected. Defaults to 500.
        AF (bool): Whether to create or not an input file for Alphafold. Default is False.

    Returns:
        pandas.DataFrame: The cleaned dataframe.
        Saves the dataframe on data_prot.csv
    """

    global gl_seq
    global gl_affinity
    global gl_cid
    global gl_smiles

    # delete all rows that contain Nan values
    pd.set_option("future.no_silent_downcasting", True)
    df = df.dropna()
    df = df.reset_index(drop=True)
    
    # delete all rows that contain sequences longer than seq_lim
    df = df[df[gl_seq].str.len() <= seq_lim]
    df = df.reset_index(drop=True)

    smiles = []
    for value in range(len(df[gl_affinity])):
        try:
            df.loc[value, gl_affinity] = float(df.loc[value, gl_affinity])
        except ValueError:
            df = df.drop([value])
    df = df.reset_index(drop=True)
    save_smiles = df[[gl_smiles, gl_cid]]
    df = (
        df.groupby([gl_seq, gl_cid,'Target Name'])[gl_affinity]
        .agg([(f"{gl_affinity} mean", "mean"), (f"{gl_affinity} sem", "sem")])
        .reset_index()
    )
    df = df.fillna(0)

    df = df[df.duplicated(subset=[gl_cid], keep=False)]
    df = df.reset_index(drop=True)
    
    for i, cid in enumerate(df[gl_cid]):
        try:
            s = save_smiles[save_smiles[gl_cid] == cid]
            s = s[gl_smiles].values[0]
            mol = Chem.MolFromSmiles(s)
            num_atoms = mol.GetNumAtoms()

            # delete rows that contain small ligands
            if num_atoms <= min_atoms or num_atoms >= max_atoms:
                df = df.drop([i])
            else:
                smiles.append(s)

        except AttributeError:
            df = df.drop([i])
    
    df = df.reset_index(drop=True)

    df["SMILES"] = smiles

    # delete values with affinity = 0
    df = df[df[f"{gl_affinity} mean"] != 0]
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)

    df["ID"] = (
        df.index
    )  # set protein IDs as the index of the dataframe as a provisional measure.

    unique_seq = np.unique(df[gl_seq])  # check the number of different proteins

    # rename all indexes if necessary

    if len(unique_seq) != len(df["ID"]):
        for i, seq in enumerate(unique_seq):
            df.loc[df[gl_seq] == seq, "ID"] = i


    return df



def mutations_from_db(db_path):
    global gl_seq
    global gl_affinity
    global gl_cid
    global gl_smiles

    df = read_db(db_path, [gl_seq, gl_affinity, gl_cid, gl_smiles,'Target Name'], display_col_names=False)

    df = clean_df(df, seq_lim=800)

    seq_similarity(df)






if __name__ == '__main__':
    gl_seq = 'BindingDB Target Chain Sequence'
    gl_smiles = 'Ligand SMILES'
    gl_affinity = 'Kd (nM)'
    gl_cid = 'PubChem CID of Ligand'

    mutations_from_db('Database_example.tsv')