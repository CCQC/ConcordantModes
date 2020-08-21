import numpy as np
import os
import shutil
import re

# Define energy search regex
molproEnergy = re.compile(r"\!RHF\-UCCSD\(T\) energy\s+(\-\d+\.\d+)")
molproSuccess = re.compile(r"Variable memory released")
psi4Energy = re.compile(r"CCSD\(T\)\s*total\s*energy\s+=\s+(\-\d+\.\d+)")
psi4Success = re.compile(r"\*\*\* Psi4 exiting successfully. Buy a developer a beer!")

rootdir = os.getcwd()
os.chdir('../zmatFiles')
with open("atomList",'r') as file:
    atoms = file.read()

os.chdir(rootdir)

atoms = atoms.split('\n')
n_atoms = len(atoms)
dispCart = "dispcart" 

os.chdir('..')
# Read in the program
if os.path.exists('prog.dat'):
    with open("prog.dat","r") as file:
        data = file.read()
    p = data.split('\n')
    progname = p[1]
else:
    progname = 'molpro'
os.chdir(rootdir)

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
    print('Prog at reap is: ')
    print(progname)
    if progname == 'molpro':
        if not re.search(molproSuccess,data):
            print('Energy failed at ' + str(i+1))
            raise RuntimeError
        energy = re.findall(molproEnergy,data)
    elif progname == 'psi4':
        if not re.search(psi4Success,data):
            print('Energy failed at ' + str(i+1))
            raise RuntimeError
        energy = re.findall(psi4Energy,data)
    else:
        print('Specified program not supported: ' + progname)
        raise RuntimeError
    Energies = np.append(Energies,energy[0])
    os.chdir('..')

outputString = ''
for i in Energies:
    outputString += i + '\n'

with open('e.dat','w') as file:
    file.write(outputString)

# with open('e.dat','r') as file:
    # energyData = file.readlines()

# print(energyData)
# for i in energyData:
    # print(i[:-2])
