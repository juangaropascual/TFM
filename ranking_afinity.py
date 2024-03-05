#Medir que ligandos estan en cada proteina
#Score esos ligandos 
#La proteína con más afinidad es...
 
import mdtraj as md
import os
from itertools import product
import numpy as np
import re
import sys

#Cargar proteina y ligandos
#prot = md.load('merge.pdb')
prot = md.load(sys. argv[1])
#path= "/home/ramon/juan/diffdock_test/results/prova4.0/ligands"
path = sys.argv[2]
distancias_A = []
distancias_B = []
i=0
def extract_number(filename):
    match = re.search(r'\d+', filename)
    return int(match.group()) if match else 0

files = os.listdir(path)
sorted_files = sorted(files, key=extract_number)
for file in sorted_files:
    
    ligand = md.load(path + '/' + file) #meter en un bucle para que se vayan analizando todos los ligandos

    top_prot, xyz_prot = prot.top, prot.xyz
    top_lig, xyz_lig = ligand.top, ligand.xyz
    top_all = top_prot.join(top_lig)
    xyz_all = np.concatenate((xyz_prot, xyz_lig), axis=1)

    prot_lig = md.Trajectory(xyz=xyz_all, topology=top_all) 
    #prot_lig.save(file + '_' + str(i) + '.pdb')  #Guarda pdb con el ligando en la proteina
    
    #Seleccionar los átomos
    chainA = prot_lig.top.select('chainid 0 name CA')
    chainB = prot_lig.top.select('chainid 1 name CA')
    ligand = prot_lig.top.select('resname UNL')

    #Calcular la distancia entre los  átomos seleccionados
    product(chainA, ligand)
    product(chainB, ligand)
    pairs_A = list(product(chainA, ligand))
    pairs_B = list(product(chainB, ligand))

    dist_A = md.compute_distances(prot_lig, pairs_A).min()
    dist_B = md.compute_distances(prot_lig, pairs_B).min()
    
    #meter en orden en una lista en que cadena está el ligando
    if dist_A < 0.4:
        distancias_A.append((i,dist_A))

    if dist_B < 0.4:
        distancias_B.append((i,dist_B))
    i+=1


for tupla in distancias_A:
    step, distance = tupla  # Desempaquetar la tupla
    if step == 0:
        print('Rank 1 ligand is in chain A')
if len(distancias_A) > len(distancias_B):
    print('Chain A has more ligands')

for tupla in distancias_B:
    step, distance = tupla  # Desempaquetar la tupla
    if step == 0:
        print('Rank 1 ligand is in chain B')
if len(distancias_B) > len(distancias_A):
    print('Chain B has more ligands')


