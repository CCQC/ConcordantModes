import numpy as np
import os
import shutil
import re

class Reap(object):
    def __init__(self, progname, zmat, dispcart, options, n_coord,eigs,indices):
        self.progname = progname
        self.zmat = zmat
        self.dispcart = dispcart
        self.options = options
        self.n_coord = n_coord
        self.dispSym = [0]
        #nate
        self.eigs = eigs
        self.indices = indices
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
        print("If something looks wrong with the final frequencies, check \
                these energies!")
        print("(Job number 1 == Reference energy) :D")

        os.chdir("./"+str(1))
        with open("output.dat",'r') as file:
            data = file.read()
        if not re.search(successRegex,data):
            print('Energy failed at ' + str('ref'))
            raise RuntimeError
        os.chdir('..')
        ref_en = float(re.findall(energyRegex,data)[0]) 
        print("Reference energy: " + str(ref_en))
        
    
        indices = self.indices
        eigs = self.eigs 
        p_en_array = np.zeros((len(eigs),len(eigs)))
        m_en_array = np.zeros((len(eigs),len(eigs)))
        rel_en_p = np.zeros((len(eigs),len(eigs)))
        rel_en_m = np.zeros((len(eigs),len(eigs)))
        relative_energies = []
       
        
        direc = 2            
        for index in indices:
            i,j = index[0], index[1] 
            p_en_array[i,j] = energy = self.reap_energies(direc,successRegex,energyRegex)
            print(energy)
            rel = ref_en - energy
            rel_en_p[i,j] = rel 
            relative_energies.append([(i,j), rel, direc]) 
            m_en_array[i,j] = energy = self.reap_energies(direc+1,successRegex,energyRegex)
            print(energy)
            rel = ref_en - energy
            rel_en_m[i,j] = rel 
            relative_energies.append([(i,j), rel, direc]) 
            direc += 2 
        
                    
        self.p_en_array = p_en_array 
        self.m_en_array = m_en_array 
        self.ref_en = ref_en
        np.set_printoptions(precision=2, linewidth = 120)
        print("Relative energies plus-displacements on the diagonal and plus/plus-displacements on the off-diagonal elements")
        print(rel_en_p)
        print("Relative energies minus-displacements on the diagonal and minus/minus-displacements on the off-diagonal elements")
        print(rel_en_m)
        os.chdir('..')
        if self.options.printout_rel_e: 
            auxillary = ''
            header = 'Index, relative energy, directory \n'
            with open("auxillary", 'a') as file:
                file.writelines(header)
            for energy in relative_energies:
                with open("auxillary", 'a') as file:
                    file.writelines(str(energy) + '\n')
    
    def reap_energies(self,direc,successRegex,energyRegex):
        os.chdir("./"+str(direc))
        with open("output.dat",'r') as file:
            data = file.read()
        if not re.search(successRegex,data):
            print('Energy failed at ' + os.getcwd())
            raise RuntimeError
        energy = float(re.findall(energyRegex,data)[0])
         
        os.chdir('..')
        return energy
        #with open("auxillary.dat", 'w') as file:
        #    for energy in relative_energies:
        #        file.writelines(energy)  
        #print(relative_energies)
        #np.set_printoptions(precision=8)
#        self.energiesDict['ref'] = float(self.Energies[0])
#        Sum = 0
#        for i in range(n_disp - np.sum(self.dispSym) - 1):
#            j = Sum % 2
#            k = Sum // 2
#            if j == 0:
#                self.energiesDict[str(k+1)+'_plus'] = float(self.Energies[i+1])
#            if j == 1:
#                self.energiesDict[str(k+1)+'_minus'] = float(self.Energies[i+1])
#            # if self.dispSym[k] == 0:
#                # if j == 0:
#                    # self.energiesDict[str(k+1)+'_plus'] = float(self.Energies[i+1])
#                # if j == 1:
#                    # self.energiesDict[str(k+1)+'_minus'] = float(self.Energies[i+1])
#            # elif self.dispSym[k] == 1:
#                # self.energiesDict[str(k+1)+'_plus'] = float(self.Energies[i+1])
#                # self.energiesDict[str(k+1)+'_minus'] = float(self.Energies[i+1])
##            Sum += 1
#            # Sum += self.dispSym[k]

   #old diagonal-only code

#        for i in range(n_disp - np.sum(self.dispSym)):
#            os.chdir("./"+str(i+1))
#            with open("output.dat",'r') as file:
#                data = file.read()
#            if not re.search(successRegex,data):
#                print('Energy failed at ' + str(i+1))
#                raise RuntimeError
#            energy = re.findall(energyRegex,data)
#            if i == 0:
#                print("Reference energy: " + str(energy[0]))
#            else:
#                print("Job number " + "{:5d}".format(i+1) + ": " + str(energy[0]))
#            if progname == 'molpro' or progname == 'psi4' or progname == 'cfour':
#                self.Energies = np.append(self.Energies,energy[0])
#            else:
#                print('Specified program not supported: ' + progname)
#                raise RuntimeError
#            os.chdir('..')
#
#        print('Relative energies')
#        for i in range(len(self.Energies)):
#            if i == 0:
#                print("Reference energy: " + str(np.float(self.Energies[i])-np.float(self.Energies[0])))
#            else:
#                print("Job number " + "{:5d}".format(i+1) + ": " + "{:+1.8f}".format(np.float(self.Energies[i])-np.float(self.Energies[0])))
#
#        self.energiesDict = {}
