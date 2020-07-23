import os
import shutil


# Read in atoms, we only really need the number of atoms
with open("atoms.txt",'r') as file:
    atoms = file.read()
atoms = atoms.split('\n')[:-1]
n_atoms = len(atoms)

# Read in dispcart
with open('dispcart','r') as file:
    disp = file.readlines()

# Initialize some dictionaries/variables
disp_dict = {}
disp_dict['disp'] = []
disp_dict['geom'] = []
n_disp = int((len(disp))/(n_atoms+1))

# Sort the displacement geometries into the dictionary and insert pertinent commas
for i in range(n_disp):
    disp_dict['disp'].append(i + 1)
    disp_dict['geom'].append([])
    for j in range(n_atoms):
        disp[i*(n_atoms + 1) + j + 1] = disp[i*(n_atoms + 1) + j + 1][:22] + ',' + disp[i*(n_atoms + 1) + j + 1][22:]
        disp[i*(n_atoms + 1) + j + 1] = disp[i*(n_atoms + 1) + j + 1][:45] + ',' + disp[i*(n_atoms + 1) + j + 1][45:]
        disp_dict['geom'][i].append(disp[i*(n_atoms + 1) + j + 1])

with open('e.dat','r') as file:
    energies = file.readlines()

# This code properly formats the output string for the e.m file
outputString = ''
outputString += 'levels = {"CCSD(T)"};\n'
outputString += 'ec = Part[Position[levels, theory], 1, 1];\n'
outputString += 'EdataX = {\n'
for i in range(n_disp):
    outputString += '\t{\n'
    outputString += ('\t\t(* This is displacement ' + str(i + 1) + ' *)\n')
    outputString += '\t\t{\n'
    for j in range(n_atoms - 1):
        outputString += ('\t\t\t{' + disp_dict['geom'][i][j][:-2] + '},\n')
    outputString += ('\t\t\t{' + disp_dict['geom'][i][n_atoms-1][:-2] + '}\n')
    outputString += '\t\t},\n'
    outputString += ('\t\t{' + energies[i][:-1] + '}\n')
    if i == n_disp - 1:
        outputString += '\t}\n'
    else:
        outputString += '\t},\n'
outputString += '};'

# Write output string to the e.m file
with open('e.m','w') as file:
    file.write(outputString)

