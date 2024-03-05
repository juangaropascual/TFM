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

top1, xyz1, time1, ul1, ua1 = p1.top, p1.xyz, p1.time, p1.unitcell_lengths, p1.unitcell_angles
top2, xyz2, time2, ul2, ua2 = p2.top, p2.xyz, p2.time, p2.unitcell_lengths, p2.unitcell_angles
top1_2 = top1.join(top2)
rotate_angles = np.linspace(0, 360, num=4, endpoint=False)  #num= nº of files with different angles

# Definir el ángulo de rotación y el eje de rotación
#angle = np.pi / 2  # 90 grados
for i in rotate_angles:
    angle = i 
    axis = np.array([1, 0, 1])  # Rotar alrededor del eje X, Z

    # Construir la matriz de rotación utilizando la fórmula de Rodrigues
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x, y, z = axis / np.linalg.norm(axis)
    rotation_matrix = np.array([[t*x*x + c,    t*x*y - z*s,  t*x*z + y*s],
                                [t*x*y + z*s,  t*y*y + c,    t*y*z - x*s],
                                [t*x*z - y*s,  t*y*z + x*s,  t*z*z + c]])

    # Aplicar la rotación a las coordenadas de la proteína
    xyz2 = np.dot(p2.xyz, rotation_matrix)

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

    p1_2 = md.Trajectory(xyz=xyz1_2, topology=top1_2, time=time1, unitcell_lengths=ul1, unitcell_angles=ua1) 
    p1_2.save(out_pdb + '_' + str(i) + '.pdb')
    print('New pdb file created: ' + out_pdb + '_' + str(i) + '.pdb')