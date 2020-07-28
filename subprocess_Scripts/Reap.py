import numpy as np
import os
import shutil
import re

# Define energy search regex
molproEnergy = re.compile(r"\!RHF\-UCCSD\(T\) energy\s+(\-\d+\.\d+)")
molproSuccess = re.compile(r"Variable memory released")

rootdir = os.getcwd()
os.chdir('../zmatFiles')
with open("atomList",'r') as file:
    atoms = file.read()

os.chdir(rootdir)

atoms = atoms.split('\n')
n_atoms = len(atoms)
dispCart = "dispcart" 

os.chdir('../mma')
with open(dispCart,'r') as file:
    disp = file.readlines()
os.chdir(rootdir)

n_disp = int((len(disp))/(n_atoms+1))
Energies = np.array([])

for i in range(n_disp):
    os.chdir("./"+str(i+1))
    with open("output.dat",'r') as file:
        data = file.read()
    if not re.search(molproSuccess,data):
        raise RuntimeError
    energy = re.findall(molproEnergy,data)
    Energies = np.append(Energies,energy[0])
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
