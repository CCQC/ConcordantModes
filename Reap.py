import numpy as np
import os
import shutil
import re

class Reap(object):
    def __init__(self, progname, zmat, dispcart, options):
        self.progname = progname
        self.zmat = zmat
        self.dispcart = dispcart
        self.options = options

    def run(self):
        # Define energy search regex
        energyRegex = re.compile(self.options.energyRegex)
        successRegex = re.compile(self.options.successRegex)

        rootdir = os.getcwd()
        n_atoms = len(self.zmat.atomList)
        
        n_disp = len(self.dispcart)
        self.Energies = np.array([])
        
        progname = self.progname
        for i in range(n_disp):
            os.chdir("./"+str(i+1))
            with open("output.dat",'r') as file:
                data = file.read()
            if not re.search(successRegex,data):
                print('Energy failed at ' + str(i+1))
                raise RuntimeError
            energy = re.findall(energyRegex,data)
            # print(i)
            # print(energy[0])
            if progname == 'molpro' or progname == 'psi4':
                self.Energies = np.append(self.Energies,energy[0])
            else:
                print('Specified program not supported: ' + progname)
                raise RuntimeError
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
