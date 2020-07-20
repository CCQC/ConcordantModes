import numpy as np
import os
import shutil
import re

# Define energy search regex
molproEnergy = re.compile(r"\!RHF\-UCCSD\(T\) energy\s+(\-\d+\.\d+)")


with open("atoms.txt",'r') as file:
    atoms = file.read()
atoms = atoms.split('\n')[:-1]

n_atoms = len(atoms)

dispCart = "dispcart" 

with open(dispCart,'r') as file:
    disp = file.readlines()

n_disp = int((len(disp))/(n_atoms+1))

os.chdir("./Disps")

Energies = np.array([])

for i in range(n_disp):
    os.chdir("./"+str(i+1))
    with open("output.dat",'r') as file:
        data = file.read()
    energy = re.findall(molproEnergy,data)
    Energies = np.append(Energies,energy[0])
    os.chdir('..')
os.chdir('..')

outputString = ''
for i in Energies:
    outputString += i + '\n'

with open('e.dat','w') as file:
    file.write(outputString)

with open('e.dat','r') as file:
    energyData = file.readlines()

print(energyData)
for i in energyData:
    print(i[:-2])
