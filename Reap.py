import numpy as np
import os
import shutil
import re

class Reap(object):
    def __init__(self, progname, zmat, dispcart):
        self.progname = progname
        self.zmat = zmat
        self.dispcart = dispcart

    def run(self):
        # Define energy search regex
        molproEnergy = re.compile(r"\!RHF\-UCCSD\(T\) energy\s+(\-\d+\.\d+)")
        molproSuccess = re.compile(r"Variable memory released")
        psi4Energy = re.compile(r"CCSD\(T\)\s*total\s*energy\s+=\s+(\-\d+\.\d+)")
        psi4Success = re.compile(r"\*\*\* Psi4 exiting successfully. Buy a developer a beer!")
        
        rootdir = os.getcwd()
        n_atoms = len(self.zmat.atomList)
        # dispCart = "dispcart" 
        
        # os.chdir(rootdir)
        
        # os.chdir('../mma')
        # with open(dispCart,'r') as file:
            # disp = file.readlines()
        # os.chdir(rootdir)
        
        n_disp = len(self.dispcart)
        self.Energies = np.array([])
        
        progname = self.progname
        for i in range(n_disp):
            os.chdir("./"+str(i+1))
            with open("output.dat",'r') as file:
                data = file.read()
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
            self.Energies = np.append(self.Energies,energy[0])
            os.chdir('..')
        
        self.energiesDict = {}
        self.energiesDict['ref'] = float(self.Energies[0])
        for i in range(n_disp - 1):
            j = i % 2
            k = i // 2
            if j == 0:
                self.energiesDict[str(k+1)+'_plus'] = float(self.Energies[i+1])
            if j == 1:
                self.energiesDict[str(k+1)+'_minus'] = float(self.Energies[i+1])

        outputString = ''
        for i in self.Energies:
            outputString += i + '\n'
        
        with open('e.dat','w') as file:
            file.write(outputString)
