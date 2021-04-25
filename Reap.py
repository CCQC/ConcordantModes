import numpy as np
import os
import shutil
import re

class Reap(object):
    def __init__(self, progname, zmat, dispcart, options, dispSym):
        self.progname = progname
        self.zmat = zmat
        self.dispcart = dispcart
        self.options = options
        self.dispSym = dispSym

    def run(self):
        """
            Define energy search regex
        """
        energyRegex = re.compile(self.options.energyRegex)
        successRegex = re.compile(self.options.successRegex)

        rootdir = os.getcwd()
        n_atoms = len(self.zmat.atomList)
        
        n_disp = len(self.dispcart)
        self.Energies = np.array([])
        
        progname = self.progname
        print("If something looks wrong with the final frequencies, check these energies!")
        print("(Job number 1 == Reference energy) :D")
        for i in range(n_disp - np.sum(self.dispSym)):
            os.chdir("./"+str(i+1))
            with open("output.dat",'r') as file:
                data = file.read()
            if not re.search(successRegex,data):
                print('Energy failed at ' + str(i+1))
                raise RuntimeError
            energy = re.findall(energyRegex,data)
            if i == 0:
                print("Reference energy: " + str(energy[0]))
            else:
                print("Job number " + "{:5d}".format(i+1) + ": " + str(energy[0]))
            if progname == 'molpro' or progname == 'psi4' or progname == 'cfour':
                self.Energies = np.append(self.Energies,energy[0])
            else:
                print('Specified program not supported: ' + progname)
                raise RuntimeError
            os.chdir('..')

        print('Relative energies')
        for i in range(len(self.Energies)):
            if i == 0:
                print("Reference energy: " + str(np.float(self.Energies[i])-np.float(self.Energies[0])))
            else:
                print("Job number " + "{:5d}".format(i+1) + ": " + "{:+1.8f}".format(np.float(self.Energies[i])-np.float(self.Energies[0])))

        self.energiesDict = {}
        self.energiesDict['ref'] = float(self.Energies[0])
        Sum = 0
        for i in range(n_disp - np.sum(self.dispSym) - 1):
            j = Sum % 2
            k = Sum // 2
            if self.dispSym[k] == 0:
                if j == 0:
                    self.energiesDict[str(k+1)+'_plus'] = float(self.Energies[i+1])
                if j == 1:
                    self.energiesDict[str(k+1)+'_minus'] = float(self.Energies[i+1])
            elif self.dispSym[k] == 1:
                self.energiesDict[str(k+1)+'_plus'] = float(self.Energies[i+1])
                self.energiesDict[str(k+1)+'_minus'] = float(self.Energies[i+1])
            Sum += 1
            Sum += self.dispSym[k]
