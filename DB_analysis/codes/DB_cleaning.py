import traceback
import numpy as np
from rdkit import Chem
import pandas as pd
from scipy.stats import sem
from difflib import SequenceMatcher
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
# uncomment the following line for interactive window (pip install ipympl):
# %matplotlib widget 
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "--display_col_names",
    type=bool,
    help="Whether or not to display the names of the columns to explore a new database. Default is False.",
    required=False,
    default=False
)
parser.add_argument(
    "--sequence",
    type=str,
    help="Name of the column containing the sequence of the proteins of the database.",
    required=False,
    default='BindingDB Target Chain Sequence'
)
parser.add_argument(
    "--threshold",
    type=float,
    help="Threshold from which the data has been studied. Default is 0.1",
    default=0.1,
    required=False,
)
parser.add_argument(
    "--affinity",
    type=str,
    help="Name of the column containing the affinity protein+ligand (can be Kd or any affinity that has the same relation with the affinity that kd: the lower the better).",
    required=False,
    default='Kd (nM)'
)
parser.add_argument(
    "--CID",
    type=str,
    help="Name of the column containing the PubChem CID (or equivalent ID) of the ligands.",
    required=False,
    default='PubChem CID'
)
parser.add_argument(
    "--SMILES",
    type=str,
    help="Name of the column containing the canonic SMILES of the ligands. ",
    required=False,
    default='Ligand SMILES'
)
parser.add_argument(
    "--match_file",
    type=str,
    help="File containing the matches with the training data of DiffDock. Default is pdbbind_match_example.csv",
    required=False,
    default='pdbbind_match_example.csv'
)
parser.add_argument(
    "--DB_file",
    type=str,
    help="Path to the database (csv or tsv) of proteins and ligands. Default is Database_example.tsv",
    required=False,
    default='Database_example.tsv'
)
args = parser.parse_args()

def read_file(file_path,columns,display_col_names=False):
    """
    Reads a tab-separated values file and returns a pandas dataframe.

    Args:
        file_path (str): The path to the file to be read.
        columns (List[str]): A list of the columns from the file that the user wants to extract.
        display_col_names (bool): Whether to display the column names. Default
        
    Returns:
        pandas.DataFrame: A pandas dataframe containing the selected data from the file.
    """
    if file_path.endswith('.tsv'):
        df = pd.read_csv(file_path, sep='\t',low_memory=False, on_bad_lines='skip')  # Use '\t' as the delimiter
    elif file_path.endswith('.csv'):
        df = pd.read_csv(file_path,low_memory=False, on_bad_lines='skip')
    else:
        raise ValueError('File format not supported. Data base must be a tsv or csv')
    
    if display_col_names:
        for col in df.columns:
            if 'Unnamed' not in col:
                print(col)
        return
    df = df[columns]
    return df

def clean_df(df, n_atoms=8, seq_lim=500, AF=False,out_file='data_prot.csv', af_out='AF_input.csv'):
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
    pd.set_option('future.no_silent_downcasting', True)
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

    save_smiles = df[[gl_smiles,gl_cid]]
    df = df.groupby([gl_seq,gl_cid])[gl_affinity].agg([(f'{gl_affinity} mean','mean'),(f'{gl_affinity} sem','sem')]).reset_index()
    df = df.fillna(0)

    df = df[df.duplicated(subset=[gl_cid],keep=False)]
    df = df.reset_index(drop=True)

    for i,cid in enumerate(df[gl_cid]):
        try:
            s = save_smiles[save_smiles[gl_cid]==cid]
            s = s[gl_smiles].values[0]
            mol = Chem.MolFromSmiles(s)
            num_atoms = mol.GetNumAtoms()

            # delete rows that contain small ligands
            if num_atoms <= n_atoms:
                df = df.drop([i])
            else:
                smiles.append(s)

        except AttributeError:
            df = df.drop([i]) 
    
    df = df.reset_index(drop=True)
    
    df['SMILES'] = smiles

    # delete values with affinity = 0
    df = df[df[f'{gl_affinity} mean']!= 0]
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)

    df['ID'] = df.index # set protein IDs as the index of the dataframe as a provisional measure.

    unique_seq = np.unique(df[gl_seq]) # check the number of different proteins

    # rename all indexes if necessary
    af = pd.DataFrame({'id':[], 'sequence':[]})
    if len(unique_seq) != len(df['ID']):
        for i,seq in enumerate(unique_seq):
            df.loc[df[gl_seq] == seq, 'ID'] = i
            
            if AF:
                row = {'id':f"'{str(i)}'", 'sequence':seq}
                af.loc[0 if pd.isnull(af.index.max()) else af.index.max() + 1] = row

 
    df.to_csv(f'output/{out_file}')

    if AF:
        af.to_csv(f'output/{af_out}')
        print(f'* {af_out} generated from {out_file} data.')

    print('')
    print(f'* Clean data saved in {out_file}')

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
    global gl_cid
    global gl_smiles

    
    un_lig = np.unique(list(df['SMILES']))
    
    for lig in range(len(un_lig)):
        sel = df[df['SMILES'] == un_lig[lig]]
        seq = list(sel[gl_seq])
        a=0
        for s1 in range(len(seq)):
            for s2 in range(s1+1, len(seq)):

                r = SequenceMatcher(None, seq[s1], seq[s2]).ratio()
                if r > 0.8:
                    print(r,s1,s2)
                    a=1
        if a==1:
            print(un_lig[lig])
            print('')

def sort_AB(df, threshold=0.1, out_A='data_A.csv', out_B='data_B.csv'):

    """
    This function sorts the dataframe by the affinity of the proteins with the same ligand.

    Args:
        df (pd.DataFrame): The dataframe to be sorted.
        threshold (float, optional): The threshold for sorting. Defaults to 0.1.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two dataframes, the rows of this dataframes can
        be paired by index to obtain two proteins that have different affinities with the same ligand. The data
        is sorted such as the affinity of the ligan with prot A is always higher than the one for B.
    """

    global gl_seq
    global gl_affinity
    global gl_cid
    global gl_smiles


    # create the new dataframes where proteins will be sorted by their affinity with the same ligand
    columns = ['Prot ID','Sequence','SMILES','Kd (nM)','kd SEM','PubChem CID']
    df_out_a = pd.DataFrame(columns=columns)
    df_out_b = pd.DataFrame(columns=columns)


    unique = np.unique(df['SMILES']) # unique ligands

    for lig in unique:

        slice = df[df['SMILES'] == lig] # obtain the slice of the dataframe that contains the ligand
        ratios = np.array(slice[f'{gl_affinity} mean'])[:, None] / np.array(slice[f'{gl_affinity} mean']) # kd_j/kd_k

        indices = np.indices(ratios.shape)
        id_slice = list(slice['ID'])

        prot = list(slice[gl_seq])
        kd = list(slice[f'{gl_affinity} mean'])
        kd_sem = list(slice[f'{gl_affinity} sem'])
        pub_cid = list(slice[gl_cid])


        mask = ratios < threshold # select only the proteins that have a ratio lower than threshold

        # since the threshold is less than 1, the least affine prot will be the one dividing
        indices_A = indices[0][mask]
        indices_B = indices[1][mask]
        ratios = ratios[mask]
        
        rows_a = []
        rows_b = []

        # Fill the new dataframes
        for i in range(len(indices_A)):

            # if len(df_out_a) == 0:
            #     df_out_a[]

            rows_a.append({'Prot ID':id_slice[indices_A[i]],
                            'Sequence':prot[indices_A[i]],
                            'SMILES':lig,
                            'Kd (nM)':kd[indices_A[i]],
                            'kd SEM': kd_sem[indices_A[i]],
                            'PubChem CID':pub_cid[indices_A[i]]})

            rows_b.append({'Prot ID':id_slice[indices_B[i]],
                            'Sequence':prot[indices_B[i]],
                            'SMILES':lig,
                            'Kd (nM)':kd[indices_B[i]],
                            'kd SEM': kd_sem[indices_B[i]],
                            'PubChem CID':pub_cid[indices_B[i]]})

        r_a = pd.DataFrame(rows_a)
        r_b = pd.DataFrame(rows_b)

        df_out_a = pd.concat([df_out_a if not df_out_a.empty else None, r_a], ignore_index=True)
        df_out_b = pd.concat([df_out_b if not df_out_b.empty else None, r_b], ignore_index=True)

    df_out_a.to_csv(f'output/{out_A}')
    df_out_b.to_csv(f'output/{out_B}')
    print(f'*** {out_A} and {out_B} successfully generated.')
    
    return df_out_a,df_out_b

def drop_if_in_training(filename,df, AF=True, out_file='clean_data.csv',af_out='AF_input.csv'):
    """
    This function removes the proteins that are included in the training of Diffusion Distance from the dataframe.

    Args:
        filename (str): The path to the file containing the proteins that were in the training data.
        df (pd.DataFrame): The dataframe containing the protein sequences and their affinities.
        AF (bool): Whether to create or not an input file for Alphafold. Default is True.

    Returns:
        pd.DataFrame: The dataframe without the proteins included in the training.

    """
    global gl_seq
    global gl_affinity
    global gl_cid
    global gl_smiles
    
    df_drop = pd.read_csv(filename)
    df = df.drop(list(df_drop['Match']))
    df = df.reset_index(drop=True)

    # If you need an AlphaFold input file:
    if AF:
        unique_seq = np.unique(df[gl_seq])
        af = pd.DataFrame({'id':[], 'sequence':[]})
        for i,seq in enumerate(unique_seq):
            id = df[df[gl_seq]==seq]['ID'].values[0]
            row = {'id':f"'{str(id)}'", 'sequence':seq}
            af.loc[0 if pd.isnull(af.index.max()) else af.index.max() + 1] = row
        
        af.to_csv(f'output/{af_out}')
        print(f'** {af_out} generated from {out_file} data.')

    df.to_csv(f'output/{out_file}')

    print(f'** Removed proteins included in training of DD. Data is in {out_file}')

    return df

if __name__ == '__main__':
    

    gl_smiles = args.SMILES
    gl_affinity = args.affinity
    gl_seq = args.sequence
    gl_cid = args.CID

    # load the database
    filename = args.DB_file
    columns = [gl_smiles, gl_affinity, gl_seq, gl_cid]
    thresh = args.threshold

    if args.display_col_names:
        df = read_file(filename,columns,display_col_names=True)
    else:   
        df = read_file(filename,columns)

        #select the useful rows of the database
        df = clean_df(df,n_atoms=8) # --> data_prot.csv can be used for searching in the training database (PDBBind) the matching sequences.

        # Print the similarity of the chains among each other for every ligand
        # seq_similarity(df) 
        file_drop = args.match_file
        try:
            # check if the selected target proteins are in the DD training databases
            df = drop_if_in_training(file_drop,df) # --> clean_data.csv can be used for creating the protein PDB files (i.e. AlphaFold)
        except:
            print(f'Skipping drop if in training since {file_drop} is not defined or is from a different dataset')
        # Sort the values of the final database in a way that kdA < kdB
        dfa,dfb = sort_AB(df,threshold=thresh) # --> data_A.csv, data_B.csv files can be paired to run the DiffDock simulation
        # print(dfa.columns)
