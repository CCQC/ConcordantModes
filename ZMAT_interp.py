import numpy as np
import os
import shutil
import re
from . import masses

class ZMAT_interp(object):
    def __init__(self):
        pass

    def run(self):
        # Define some regexes
        zmatBeginRegex  = re.compile(r"ZMAT begin")
        zmatEndRegex    = re.compile(r"ZMAT end")
        firstAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s+\n")
        secondAtomRegex = re.compile("^\s*([A-Z][a-z]*)\s(\d+)\s([A-Za-z0-9_]+)\s*\n")
        thirdAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s(\d+)\s([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        fullAtomRegex   = re.compile("^\s*([A-Z][a-z]*)\s(\d+)\s([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        variableRegex   = re.compile("^\s*([A-Za-z0-9_]+)\s*\=\s*(-?\d+\.\d+)\s*\n")
        cartBeginRegex  = re.compile(r"cart begin")
        cartEndRegex    = re.compile(r"cart end")
        # cartesianRegex  = re.compile("^\s+[A-Z][A-Za-z]*\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\n")
        cartesianRegex  = re.compile("[A-Z][A-Za-z]*\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\n")
        
        # Read in the ZMAT file
        with open("zmat",'r') as file:
            output = file.readlines()
        
        tempOutput = output.copy()
        
        # Slice out the ZMAT from the input
        zmatRange = []
        
        for i in range(len(output)):
            begZMAT = re.search(zmatBeginRegex,output[i])
            if begZMAT:
                zmatRange.append(i)
        
        for i in range(len(output) - zmatRange[0]):
            endZMAT = re.search(zmatEndRegex,output[i + zmatRange[0]])
            if endZMAT:
                zmatRange.append(i + zmatRange[0])
                break
        
        output = output[zmatRange[0]:zmatRange[1]]
        
        # Initialize necessary lists
        self.atomList           = []
        self.bondIndices        = []
        self.bondVariables      = []
        self.angleIndices       = []
        self.angleVariables     = []
        self.torsionIndices     = []
        self.torsionVariables   = []
        self.variableDictionary = {}
        # Now to reap the ZMAT data 
        for i in range(len(output)):
            # This case if we are at the first atom of the ZMAT
            if re.search(firstAtomRegex,output[i]):
                self.atomList.append(re.findall(firstAtomRegex,output[i])[0])
                firstIndex = i
            # Second atom of the ZMAT, will have one bond term
            if re.search(secondAtomRegex,output[i]):
                List = re.findall(secondAtomRegex,output[i])[0]
                self.atomList.append(List[0])
                self.bondIndices.append([str(i-firstIndex+1),List[1]])
                self.bondVariables.append(List[2])
            # Third atom of the ZMAT, will have bond and angle term
            if re.search(thirdAtomRegex,output[i]):
                List = re.findall(thirdAtomRegex,output[i])[0]
                self.atomList.append(List[0])
                self.bondIndices.append([str(i-firstIndex+1),List[1]])
                self.bondVariables.append(List[2])
                self.angleIndices.append([str(i-firstIndex+1),List[1],List[3]])
                self.angleVariables.append(List[4])
            # All remaining ZMAT atoms, will have bond, angle, and torsion term
            if re.search(fullAtomRegex,output[i]):
                List = re.findall(fullAtomRegex,output[i])[0]
                self.atomList.append(List[0])
                self.bondIndices.append([str(i-firstIndex+1),List[1]])
                self.bondVariables.append(List[2])
                self.angleIndices.append([str(i-firstIndex+1),List[1],List[3]])
                self.angleVariables.append(List[4])
                self.torsionIndices.append([str(i-firstIndex+1),List[1],List[3],List[5]])
                self.torsionVariables.append(List[6])
            # Now to reap the variables and their values
            if re.search(variableRegex,output[i]):
                List = re.findall(variableRegex,output[i])[0]
                self.variableDictionary[List[0]] = float(List[1])

        # print('Final ZMAT Info:')
        # print(self.atomList)
        # print(self.bondIndices)
        # print(self.bondVariables)
        # print(self.angleIndices)
        # print(self.angleVariables)
        # print(self.torsionIndices)
        # print(self.torsionVariables)
        # print(self.variableDictionary)
        """ This code may be helpful if I make my own ZMAT generator. """
        
        cartRange = []
        
        for i in range(len(tempOutput)):
            begCart = re.search(cartBeginRegex,tempOutput[i])
            if begCart:
                cartRange.append(i)
                break
        
        for i in range(len(tempOutput) - cartRange[0]):
            endCart = re.search(cartEndRegex,tempOutput[i + cartRange[0]])
            if endCart:
                cartRange.append(i + cartRange[0])
                break

        tempOutput = tempOutput[cartRange[0]:cartRange[1]] 
        # Find the cartesian coordinates, in bohr
        self.Cartesians = []
        for i in range(len(tempOutput)):
            if re.search(cartesianRegex,tempOutput[i]):
                temp = re.findall(cartesianRegex,tempOutput[i])
                self.Cartesians.append(temp[0])
        self.Cartesians = np.array(self.Cartesians).astype(float)
        # print(self.Cartesians)
        
        self.masses = [masses.get_mass(label) for label in self.atomList]
