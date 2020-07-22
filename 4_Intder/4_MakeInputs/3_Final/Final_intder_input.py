import numpy as np
import os
import shutil
import re

# Do some directory tomfoolery here to get all the necessary files
with open('intder_Final_template.dat','r') as file:
    template = file.readlines()

root = os.getcwd()
if os.path.exists(root + '/zmatFiles'):
    shutil.rmtree(root + '/zmatFiles')

os.chdir(root + '/..')
shutil.copytree('4_zmatFiles',root + '/zmatFiles')
os.chdir(root + '/zmatFiles')

# Read in the ZMAT data
with open('atomList','r') as file:
    atoms = file.read()
with open('bondIndices','r') as file:
    bonds = file.read()
with open('angleIndices','r') as file:
    angles = file.read()
with open('torsionIndices','r') as file:
    torsions = file.read()
with open('Cartesians','r') as file:
    cartesians = file.read()

os.chdir(root)

# split the data into arrays

# print(atoms.split('\n'))
atoms = atoms.split('\n') 
# print(bonds.split('\n'))
bonds = bonds.split('\n')
for i in range(len(bonds)):
    bonds[i] = bonds[i].split(' ')
# print(angles.split('\n'))
angles = angles.split('\n')
for i in range(len(angles)):
    angles[i] = angles[i].split(' ')
# print(torsions.split('\n'))
torsions = torsions.split('\n') 
for i in range(len(torsions)):
    torsions[i] = torsions[i].split(' ')
# print(cartesians.split('\n'))
cartesians = cartesians.split('\n')
for i in range(len(cartesians)):
    cartesians[i] = cartesians[i].split(' ')

# Now we can build the intder.inp file
template[-1] = '{:5d}'.format(1) + template[-1][35:]
template[-1] = '{:5d}'.format(2000) + template[-1]
template[-1] = '{:5d}'.format(0) + template[-1]
template[-1] = '{:5d}'.format(2) + template[-1]
template[-1] = '{:5d}'.format(0) + template[-1]
template[-1] = '{:5d}'.format(len(atoms)*3 - 6) + template[-1]
template[-1] = '{:5d}'.format(len(atoms)) + template[-1]

# Write the bond internal coords
if bonds[0][0] != '':
    for i in bonds:
        buffString = ''
        buffString += 'STRE '
        for j in i:
            buffString += '{:5d}'.format(int(j))
        buffString += '\n'
        template.append(buffString)

# Write the angle internal coords
if angles[0][0] != '':
    for i in angles:
        buffString = ''
        buffString += 'BEND '
        for j in i:
            buffString += '{:5d}'.format(int(j))
        buffString += '\n'
        template.append(buffString)

# Write the torsion internal coords
if torsions[0][0] != '':
    for i in torsions:
        buffString = ''
        buffString += 'TORS '
        for j in i:
            buffString += '{:5d}'.format(int(j))
        buffString += '\n'
        template.append(buffString)

# Now we can write in the cartesian coordinates
if len(cartesians) > 0:
    for i in cartesians:
        buffString = ''
        buffString += '{:19.10f}'.format(float(i[0]))
        buffString += '{:21.10f}'.format(float(i[1]))
        buffString += '{:20.10f}'.format(float(i[2]))
        buffString += '\n'
        template.append(buffString)

# And finally, print out the atoms
if len(atoms) > 0:
    buffString = ''
    for i in atoms:
        buffString += '{:>12s}'.format(i)
        if len(buffString) > 70:
            buffString += '\n'
            template.append(buffString)
            buffString = ''
    if len(buffString) > 0:
        # buffString += '\n'
        template.append(buffString)
    # This statement ensures there is no extra line break at the end of the file if the atom
    # number is a multiple of 6.
    else:
        template[-1] = template[-1][:-1]

with open('intder.inp','w') as file:
    file.writelines(template)



