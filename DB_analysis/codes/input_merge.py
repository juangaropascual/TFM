import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "--chain_A",
    type=str,
    help="CSV file with the Chain A information",
    required=True,
)
parser.add_argument(
    "--chain_B",
    type=str,
    help="CSV file with the Chain B information",
    required=True,
)
parser.add_argument(
    "--name",
    type=str,
    help="Name of the run",
    required=True,
)
args = parser.parse_args()


# Archivos de entrada y salida
data_a = pd.read_csv(args.chain_A)
data_b = pd.read_csv(args.chain_B)
run_name = args.name
input_DD = {'complex_name':[], 'chainA': [], 'chainB': [], 'SMILES': []}
path = 'output/prot_structures/' + run_name + '/'

complex_name = data_a['Prot ID'].astype(str) + '_' + data_b['Prot ID'].astype(str) + '_' + data_a['PubChem CID'].astype(str)

input_DD['complex_name'] = complex_name
input_DD['chainA'] = path + data_a['Prot ID'].astype(str) + '.pdb'
input_DD['chainB'] = path + data_b['Prot ID'].astype(str) + '.pdb'
input_DD['SMILES'] = data_a['SMILES']

df = pd.DataFrame(input_DD)
df.to_csv('output/merge_input.csv', index=False)

