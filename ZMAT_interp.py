import numpy as np
import os
import shutil
import re
from . import masses

# class ZMAT_interp(object):
class ZMAT(object):
    def __init__(self):
        self.amu_elMass = 5.48579909065*(10**(-4))

    def run(self):
        # Define some regexes
        zmatBeginRegex  = re.compile(r"ZMAT begin")
        zmatEndRegex    = re.compile(r"ZMAT end")
        firstAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s*\n")
        secondAtomRegex = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        thirdAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        fullAtomRegex   = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        variableRegex   = re.compile("\s*([A-Za-z0-9_]+)\s*\=\s*(-?\d+\.\d+)\s*\n")
        cartBeginRegex  = re.compile(r"cart begin")
        cartEndRegex    = re.compile(r"cart end")
        cartesianRegex  = re.compile("[A-Z][A-Za-z]*\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\n")
        dividerRegex    = re.compile("^\s*\-\-\-\s*\n")
        
        # Read in the ZMAT file
        with open("zmat",'r') as file:
            output = file.readlines()
        
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
        
        zmatOutput = output[zmatRange[0]:zmatRange[1]].copy()
        
        for i in range(len(zmatOutput)):
            divider = re.search(dividerRegex,zmatOutput[i])
            if divider:
                divideIndex = i
                break
        if divider:
            zmatOutputInit  = zmatOutput[:divideIndex].copy()
            zmatOutputFinal = zmatOutput[divideIndex+1:].copy()
        else:
            zmatOutputInit  = zmatOutput.copy()

        # Initialize necessary lists
        self.atomList                = []
        self.bondIndices             = []
        self.bondVariables           = []
        self.angleIndices            = []
        self.angleVariables          = []
        self.torsionIndices          = []
        self.torsionVariables        = []
        self.variableDictionaryInit  = {}
        self.variableDictionaryFinal = {}

        # Now to reap the ZMAT data 
        count = 0
        for i in range(len(zmatOutputInit)):
            # This case if we are at the first atom of the ZMAT
            if re.search(firstAtomRegex,zmatOutputInit[i]) and count < 1:
                self.atomList.append(re.findall(firstAtomRegex,zmatOutputInit[i])[0])
                firstIndex = i
                count += 1
            # Second atom of the ZMAT, will have one bond term
            if re.search(secondAtomRegex,zmatOutputInit[i]):
                List = re.findall(secondAtomRegex,zmatOutputInit[i])[0]
                self.atomList.append(List[0])
                self.bondIndices.append([str(i-firstIndex+1),List[1]])
                self.bondVariables.append(List[2])
            # Third atom of the ZMAT, will have bond and angle term
            if re.search(thirdAtomRegex,zmatOutputInit[i]):
                List = re.findall(thirdAtomRegex,zmatOutputInit[i])[0]
                self.atomList.append(List[0])
                self.bondIndices.append([str(i-firstIndex+1),List[1]])
                self.bondVariables.append(List[2])
                self.angleIndices.append([str(i-firstIndex+1),List[1],List[3]])
                self.angleVariables.append(List[4])
            # All remaining ZMAT atoms, will have bond, angle, and torsion term
            if re.search(fullAtomRegex,zmatOutputInit[i]):
                List = re.findall(fullAtomRegex,zmatOutputInit[i])[0]
                self.atomList.append(List[0])
                self.bondIndices.append([str(i-firstIndex+1),List[1]])
                self.bondVariables.append(List[2])
                self.angleIndices.append([str(i-firstIndex+1),List[1],List[3]])
                self.angleVariables.append(List[4])
                self.torsionIndices.append([str(i-firstIndex+1),List[1],List[3],List[5]])
                self.torsionVariables.append(List[6])
            # Now to reap the variables and their values
            if re.search(variableRegex,zmatOutputInit[i]):
                List = re.findall(variableRegex,zmatOutputInit[i])[0]
                self.variableDictionaryInit[List[0]] = float(List[1])

        
        if divider:
            for i in range(len(zmatOutputFinal)):
                if re.search(variableRegex,zmatOutputFinal[i]):
                    List = re.findall(variableRegex,zmatOutputFinal[i])[0]
                    self.variableDictionaryFinal[List[0]] = float(List[1])
        else:
            self.variableDictionaryFinal = self.variableDictionaryInit.copy()
        """ 
            And now we must temper the torsion angles! For consistency's sake we will
            force them to lie between -90 deg and +270 deg.
        """
        """ Handle Variable lists separately. First the INIT: """
        for i in range(len(self.torsionVariables)):
            Condition1 = float(self.variableDictionaryInit[self.torsionVariables[i]]) <= -90.
            Condition2 = float(self.variableDictionaryInit[self.torsionVariables[i]]) >= 270.
            buff = np.floor(abs(float(self.variableDictionaryInit[self.torsionVariables[i]]))/360)
            if Condition1:
                self.variableDictionaryInit[self.torsionVariables[i]] = float(self.variableDictionaryInit[self.torsionVariables[i]])
                self.variableDictionaryInit[self.torsionVariables[i]] += 360.*buff
                if float(self.variableDictionaryInit[self.torsionVariables[i]]) <= -90.:
                    self.variableDictionaryInit[self.torsionVariables[i]] += 360.
            if Condition2:
                self.variableDictionaryInit[self.torsionVariables[i]] = float(self.variableDictionaryInit[self.torsionVariables[i]])
                self.variableDictionaryInit[self.torsionVariables[i]] -= 360.*buff
                if float(self.variableDictionaryInit[self.torsionVariables[i]]) >= 270.:
                    self.variableDictionaryInit[self.torsionVariables[i]] -= 360.
        """ Then the Final. This can probably be structured more elegantly, but this works and isn't too computationally demanding. """
        for i in range(len(self.torsionVariables)):
            Condition1 = float(self.variableDictionaryFinal[self.torsionVariables[i]]) <= -90.
            Condition2 = float(self.variableDictionaryFinal[self.torsionVariables[i]]) >= 270.
            buff = np.floor(abs(float(self.variableDictionaryFinal[self.torsionVariables[i]]))/360)
            if Condition1:
                self.variableDictionaryFinal[self.torsionVariables[i]] = float(self.variableDictionaryFinal[self.torsionVariables[i]])
                self.variableDictionaryFinal[self.torsionVariables[i]] += 360.*buff
                if float(self.variableDictionaryFinal[self.torsionVariables[i]]) <= -90.:
                    self.variableDictionaryFinal[self.torsionVariables[i]] += 360.
            if Condition2:
                self.variableDictionaryFinal[self.torsionVariables[i]] = float(self.variableDictionaryFinal[self.torsionVariables[i]])
                self.variableDictionaryFinal[self.torsionVariables[i]] -= 360.*buff
                if float(self.variableDictionaryFinal[self.torsionVariables[i]]) >= 270.:
                    self.variableDictionaryFinal[self.torsionVariables[i]] -= 360.
       
        """
            This code slices out the Cartesian coordinates
        """
        cartRange = []
        
        for i in range(len(output)):
            begCart = re.search(cartBeginRegex,output[i])
            if begCart:
                cartRange.append(i)
                break
        
        for i in range(len(output) - cartRange[0]):
            endCart = re.search(cartEndRegex,output[i + cartRange[0]])
            if endCart:
                cartRange.append(i + cartRange[0])
                break

        cartOutput = output[cartRange[0]:cartRange[1]].copy()
        
        for i in range(len(cartOutput)):
            divider = re.search(dividerRegex,cartOutput[i])
            if divider:
                divideIndex = i
                break
        if divider:
            cartOutputInit  = cartOutput[:divideIndex].copy()
            cartOutputFinal = cartOutput[divideIndex+1:].copy()
        else:
            cartOutputInit  = cartOutput.copy()
        
        # Find the cartesian coordinates, in bohr
        self.CartesiansInit = []
        self.CartesiansFinal = []
        for i in range(len(cartOutputInit)):
            if re.search(cartesianRegex,cartOutputInit[i]):
                temp = re.findall(cartesianRegex,cartOutputInit[i])
                self.CartesiansInit.append(temp[0])
        self.CartesiansInit = np.array(self.CartesiansInit).astype(float)

        if divider:
            for i in range(len(cartOutputFinal)):
                if re.search(cartesianRegex,cartOutputFinal[i]):
                    temp = re.findall(cartesianRegex,cartOutputFinal[i])
                    self.CartesiansFinal.append(temp[0])
            self.CartesiansFinal = np.array(self.CartesiansFinal).astype(float)
        else:
            self.CartesiansFinal = self.CartesiansInit.copy()

        """
            The masses are assigned to the respective atom from the masses.py file
        """
        self.masses = [masses.get_mass(label) for label in self.atomList]
        for i in range(len(self.masses)):
            self.masses[i] = self.masses[i]/self.amu_elMass
