import pandas as pd

# Archivos de entrada y salida
data_a = pd.read_csv('data_A.csv')
data_b = pd.read_csv('data_B.csv')

input_DD = {'complex_name':[], 'chainA': [], 'chainB': [], 'SMILES': []}
path = '/home/ramon/juan/structures/kd_articles/AF_structures/'

complex_name = data_a['Prot ID'].astype(str) + '_' + data_b['Prot ID'].astype(str) + '_' + data_a['PubChem CID'].astype(str)

input_DD['complex_name'] = complex_name
input_DD['chainA'] = path + data_a['Prot ID'].astype(str) + '.pdb'
input_DD['chainB'] = path + data_b['Prot ID'].astype(str) + '.pdb'
input_DD['SMILES'] = data_a['SMILES']

df = pd.DataFrame(input_DD)
df.to_csv('input_DD.csv', index=False)

