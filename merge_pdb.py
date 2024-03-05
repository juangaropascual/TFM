import mdtraj as md
import numpy as np
import sys

#Recognise each argument
for i in range(0, len(sys.argv)):
    if sys. argv[i] == "--chainA":
        chainA = sys.argv[i + 1]
    elif sys.argv[i] == "--chainB":
        chainB = sys.argv[i + 1]
    elif sys.argv[i] == "--out":
        out_pdb = sys.argv[i + 1]
    if len(sys.argv) < 7:
        print('Error: python coordinates.py --chainA chainA.pdb --chainB chainB.pdb --out out_file')
        exit()

#Load the structures
p1 = md.load(chainA) 
p2 = md.load(chainB)

top1, xyz1 = p1.top, p1.xyz
top2, xyz2 = p2.top, p2.xyz
top1_2 = top1.join(top2)

#Asegurarse que las estructuras estan separadas
if np.linalg.norm(xyz1.mean(axis = 1) - xyz2.mean(axis = 1)) < 10:
    mean1 = xyz1.mean(axis = 1)
    mean2 = xyz2.mean(axis = 1)
    vector = mean2 - mean1
    distance = np.linalg.norm(vector)
    unit_vector = vector / distance
    separation_vector = unit_vector * 10
    xyz2 += separation_vector

xyz1_2 = np.concatenate((xyz1, xyz2), axis=1)

p1_2 = md.Trajectory(xyz=xyz1_2, topology=top1_2) 
p1_2.save(out_pdb + '.pdb')
print('New pdb file created: ' + out_pdb + '.pdb')