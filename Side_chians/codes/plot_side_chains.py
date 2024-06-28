import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm


fold = 'AF'
# Read data
data_1 = pd.read_csv(f'results/docking_side-chains_{fold}.csv')
data_2 = pd.read_csv(f'results/docking_side-chains-2_{fold}.csv')
data_3 = pd.read_csv(f'results/docking_side-chains-3_{fold}.csv')

# Agrupar por 'ID'
grouped_1 = data_1.groupby('ID')
grouped_2 = data_2.groupby('ID')
grouped_3 = data_3.groupby('ID')

# Crear una figura y un eje
fig, ax = plt.subplots()

# IDs deseados
desired_ids = ['Backbone-dependent (2010)',
                'Backbone-dependent (2002)',
                'Backbone-independent']
colors = {id: plt.cm.plasma(i/len(desired_ids)) for i, id in enumerate(desired_ids)}

# Crear un diccionario para almacenar las líneas de la leyenda y sus etiquetas

legend_lines = {}

# Para cada grupo, trazar la 'Confidence' frente al número de la primera columna
for name, group in grouped_1:
    if name in desired_ids:
        group = group[group['Confidence'] >= -10]  # Filtrar los datos
        line, = ax.plot(group['Unnamed: 0'], group['Confidence'], label=name, color=colors[name])
        legend_lines[name] = line

for name, group in grouped_2:
    if name in desired_ids:
        group = group[group['Confidence'] >= -10]  # Filtrar los datos
        line, = ax.plot(group['Unnamed: 0'], group['Confidence'], label=name, color=colors[name])
        legend_lines[name] = line

for name, group in grouped_3:
    if name in desired_ids:
        group = group[group['Confidence'] >= -10]  # Filtrar los datos
        line, = ax.plot(group['Unnamed: 0'], group['Confidence'], label=name, color=colors[name])
        legend_lines[name] = line

# Configurar el gráfico
ax.set_xlabel('Rank')
ax.set_ylabel('Confidence')
ax.set_title(f'Confidence comparation  ({fold})')
ax.set_ylim(-5, 1)

# Crear la leyenda a partir del diccionario de líneas
ax.legend(legend_lines.values(), legend_lines.keys())
# print(legend_lines)
# Mostrar el gráfico
plt.savefig(f'results/figures/side-chains_new_{fold}.svg')