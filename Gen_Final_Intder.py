import numpy as np
import subprocess as sp
import re
import os
import shutil


class Final_Intder(object):
    def __init__(self):
        self.rootdir = os.getcwd()
        
        # Define Regexes
        self.fcRegex = re.compile('\s+\d+\s+\d+\s+(-?\d+\.\d+)')
        self.f16Regex = re.compile('\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)')

    def run(self):
        # Open some files and read in the pertinent values
        os.chdir('../mma')
        with open('fc','r') as file:
            fcdata = file.readlines()
        os.chdir(self.rootdir)
        
        with open('file16','r') as file:
            f16 = file.readlines()
        
        shutil.copy('file16','old_file16')
        
        f16_header = f16[0]
        
        # generate list of diagonal force constants
        fc = np.array([])
        for i in fcdata:
            fc = np.append(fc,re.findall(self.fcRegex,i))
        
        # generate list of lower level diagonal + off-diagonal force constants
        f16_arr = np.array([])
        for i in f16:
            f16_arr = np.append(f16_arr,re.findall(self.f16Regex,i))
        
        dim = int(np.sqrt(len(f16_arr)))
        
        f16_mat = np.reshape(f16_arr,(dim,dim))
        
        for i in range(dim):
            # This is a good debugging check, force constants should be similar
            # print(f16_new[i][i])
            # print(fc[i])
            f16_mat[i][i] = fc[i]
        
        f16_flat = f16_mat.flatten()
        
        f16_outputString = f16_header
        for i in range(len(f16_flat)):
            if (i+1)%3 == 0:
                f16_outputString += "{:20.10f}".format(float(f16_flat[i])) + '\n'
            else:
                f16_outputString += "{:20.10f}".format(float(f16_flat[i]))
        
        f16_outputString = f16_outputString[:-1]
        
        # This code will be for generating all the necessary intder files for the last frequency computation
        with open('file16','w') as file:
            file.write(f16_outputString)
        os.system('/home/vulcan/mel64643/bin/MixedHessian/Temporary_Scripts/INTDER < intder.inp > intder.out')

