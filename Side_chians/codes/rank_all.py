import glob
import pandas as pd

run_name = 'all'
# Obtén una lista de todos los archivos que comienzan con 'docking_{run_name}_{fold}.csv'
file_list = glob.glob('/home/ramon/juan/CODI/single_prot/results/docking_*_*.csv')

# Divide los archivos en tres listas basadas en el valor de {fold}
af_files = [file for file in file_list if '_AF.' in file]
df_files = [file for file in file_list if '_DF.' in file]
of_files = [file for file in file_list if '_OF.' in file]

# Lee cada archivo en un DataFrame y guárdalos en una lista
af_dfs = [pd.read_csv(file) for file in af_files]
df_dfs = [pd.read_csv(file) for file in df_files]
of_dfs = [pd.read_csv(file) for file in of_files]

# Concatena todos los DataFrames en la lista y guarda cada uno en un archivo
af_all = pd.concat(af_dfs)
df_all = pd.concat(df_dfs)
# of_all = pd.concat(of_dfs)

af_all.to_csv('results/docking_all_AF.csv', index=False)
df_all.to_csv('results/docking_all_DF.csv', index=False)
# of_all.to_csv('results/docking_all_OF.csv', index=False)

# Concatena los tres archivos en uno
all_df = pd.concat([af_all, df_all])
all_df.to_csv('results/docking_all_all.csv', index=False)



# Obtén una lista de todos los archivos que comienzan con 'most_affine_*.csv'
file_list = glob.glob(f'/home/ramon/juan/CODI/single_prot/results/most_affine_*.csv')

# Lee cada archivo en un DataFrame, normaliza la columna 'Kd' y guárdalos en una lista
dfs = []
for file in file_list:
    df = pd.read_csv(file)
    df['Kd'] = df['Kd'] / df['Kd'].max()
    dfs.append(df)

# Concatena todos los DataFrames en la lista y guarda el resultado en un archivo
all_df = pd.concat(dfs)
print(all_df)
all_df.to_csv(f'results/most_affine_{run_name}.csv', index=False)