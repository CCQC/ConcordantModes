import numpy as np
import os
import shutil
import re
from . import masses

class ZMAT(object):
    def __init__(self):
        self.amu_elMass = 5.48579909065*(10**(-4))

    def run(self):
        # Define some regexes
        zmatBeginRegex  = re.compile(r"Input from ZMAT file")
        zmatEndRegex1   = re.compile(r"CFOUR")
        zmatEndRegex2   = re.compile(r"ACESII")
        firstAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s*\n")
        secondAtomRegex = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        thirdAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        fullAtomRegex   = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        variableRegex   = re.compile("^\s*([A-Za-z0-9_]+)\s*\=\s*(-?\d+\.\d+)\s*\n")
        cartesianRegex  = re.compile("^\s+[A-Z][A-Za-z]*\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\n")
        
        
        # Read in the output file
        with open("output.dat",'r') as file:
            output = file.readlines()
        
        tempOutput = output.copy()
        
        # Slice out the ZMAT from the output
        zmatRange = []
        
        for i in range(len(output)):
            begZMAT = re.search(zmatBeginRegex,output[i])
            if begZMAT:
                zmatRange.append(i)
        
        for i in range(len(output) - zmatRange[0]):
            endZMAT1 = re.search(zmatEndRegex1,output[i + zmatRange[0]])
            endZMAT2 = re.search(zmatEndRegex2,output[i + zmatRange[0]])
            if endZMAT1 or endZMAT2:
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
        
        """ 
            And now we must temper the torsion angles! For consistency's sake we will
            force them to lie between -90 deg and +270 deg.
        """
        for i in range(len(self.torsionVariables)):
            Condition1 = float(self.variableDictionary[self.torsionVariables[i]]) <= -90.
            Condition2 = float(self.variableDictionary[self.torsionVariables[i]]) >= 270.
            buff = np.floor(abs(float(self.variableDictionary[self.torsionVariables[i]])/360))
            if Condition1:
                self.variableDictionary[self.torsionVariables[i]] = float(self.variableDictionary[self.torsionVariables[i]])
                self.variableDictionary[self.torsionVariables[i]] += 360.*buff
                if float(self.variableDictionary[self.torsionVariables[i]]) <= -90.:
                    self.variableDictionary[self.torsionVariables[i]] += 360.
            elif Condition2:
                self.variableDictionary[self.torsionVariables[i]] = float(self.variableDictionary[self.torsionVariables[i]])
                self.variableDictionary[self.torsionVariables[i]] -= 360.*buff
                if float(self.variableDictionary[self.torsionVariables[i]]) >= 270.:
                    self.variableDictionary[self.torsionVariables[i]] -= 360.
        
        # Make a cartesian output slice
        for i in range(len(tempOutput)):
            if re.search(cartesianRegex,tempOutput[i]):
                cartIndex = i
                break
        
        # Find the cartesian coordinates, in bohr
        self.Cartesians = []
        for i in range(len(tempOutput)):
            if re.search(cartesianRegex,tempOutput[i]):
                temp = re.findall(cartesianRegex,tempOutput[i])
                self.Cartesians.append(temp[0])
            if i > (cartIndex + len(self.atomList) + 1):
                break

        """
            Convert the angle and torsion variables to np arrays and remove duplicates
        """
        flatAngles = np.array(self.angleVariables)
        flatAngles = np.unique(flatAngles)
        flatTorsions = np.array(self.torsionVariables)
        flatTorsions = np.unique(flatTorsions)

        """
            Convert angles and torsions to radians
        """
        for i in flatAngles:
            self.variableDictionary[i] = float(self.variableDictionary[i])*(np.pi/180)
        for i in flatTorsions:
            self.variableDictionary[i] = float(self.variableDictionary[i])*(np.pi/180)
       
        """
            The proper masses are pulled from the masses.py file
        """
        self.masses = [masses.get_mass(label) for label in self.atomList]
        for i in range(len(self.masses)):
            self.masses[i] = self.masses[i]/self.amu_elMass
