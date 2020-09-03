import numpy as np
import os
import shutil
import re

class ZMAT(object):
    def __init__(self):
        pass

    def run(self):
        # Define some regexes
        zmatBeginRegex  = re.compile(r"Input from ZMAT file")
        zmatEndRegex1   = re.compile(r"CFOUR")
        zmatEndRegex2   = re.compile(r"ACESII")
        firstAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s+\n")
        secondAtomRegex = re.compile("^\s*([A-Z][a-z]*)\s(\d+)\s([A-Za-z0-9_]+)\s*\n")
        thirdAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s(\d+)\s([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
        fullAtomRegex   = re.compile("^\s*([A-Z][a-z]*)\s(\d+)\s([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n")
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
                self.variableDictionary[List[0]] = List[1]
        
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
        
        atomListOutput           = ''
        bondIndicesOutput        = ''
        bondVariablesOutput      = ''
        angleIndicesOutput       = ''
        angleVariablesOutput     = ''
        torsionIndicesOutput     = ''
        torsionVariablesOutput   = ''
        variableDictionaryOutput = ''
        cartesianOutput          = ''
        
        # Write the atom list output
        if len(self.atomList) > 0:
            for i in range(len(self.atomList)-1):
                atomListOutput += self.atomList[i] + '\n'
            atomListOutput += self.atomList[len(self.atomList)-1]
        
        # Write the bond outputs
        if len(self.bondIndices) > 0:
            for i in range(len(self.bondIndices)-1):
                bondIndicesOutput += self.bondIndices[i][0] + ' ' + self.bondIndices[i][1] + '\n'
            bondIndicesOutput += self.bondIndices[len(self.bondIndices)-1][0] + ' ' + self.bondIndices[len(self.bondIndices)-1][1]
            for i in range(len(self.bondVariables)-1):
                bondVariablesOutput += self.bondVariables[i] + '\n'
            bondVariablesOutput += self.bondVariables[len(self.bondVariables)-1]
        
        # Write the angle outputs
        if len(self.angleIndices) > 0:
            for i in range(len(self.angleIndices)-1):
                angleIndicesOutput += self.angleIndices[i][0] + ' ' + self.angleIndices[i][1] + ' ' + self.angleIndices[i][2] + '\n'
            angleIndicesOutput += self.angleIndices[len(self.angleIndices)-1][0] + ' ' + self.angleIndices[len(self.angleIndices)-1][1] + ' ' + self.angleIndices[len(self.angleIndices)-1][2]
            for i in range(len(self.angleVariables)-1):
                angleVariablesOutput += self.angleVariables[i] + '\n'
            angleVariablesOutput += self.angleVariables[len(self.angleVariables)-1]
        
        # Write the torsion outputs
        if len(self.torsionIndices) > 0:
            for i in range(len(self.torsionIndices)-1):
                torsionIndicesOutput += self.torsionIndices[i][0] + ' ' + self.torsionIndices[i][1] + ' ' + self.torsionIndices[i][2] + ' ' + self.torsionIndices[i][3] + '\n'
            torsionIndicesOutput += self.torsionIndices[len(self.torsionIndices)-1][0] + ' ' + self.torsionIndices[len(self.torsionIndices)-1][1] + ' ' + self.torsionIndices[len(self.torsionIndices)-1][2] + ' ' + self.torsionIndices[len(self.torsionIndices)-1][3]
            for i in range(len(self.torsionVariables)-1):
                torsionVariablesOutput += self.torsionVariables[i] + '\n'
            torsionVariablesOutput += self.torsionVariables[len(self.torsionVariables)-1]
        # Write out the variables and their values
        # self.bondVariables = list(dict.fromkeys(self.bondVariables))
        # self.angleVariables = list(dict.fromkeys(self.angleVariables))
        # self.torsionVariables = list(dict.fromkeys(self.torsionVariables))
        
        flatAngles = np.array(self.angleVariables)
        flatTorsions = np.array(self.torsionVariables)
        for i in flatAngles:
            self.variableDictionary[i] = str(float(self.variableDictionary[i]) * (np.pi/180))
        for i in flatTorsions:
            self.variableDictionary[i] = str(float(self.variableDictionary[i]) * (np.pi/180))
        
        # variables = np.array(self.bondVariables)
        # variables = np.append(variables,self.angleVariables)
        # variables = np.append(variables,self.torsionVariables)
        # if len(self.variableDictionary) > 0:
            # if len(variables) > 0:
                # for i in range(len(variables)-1):
                    # variableDictionaryOutput += variables[i] + ' ' + self.variableDictionary[variables[i]] + '\n'
                # variableDictionaryOutput += variables[len(variables)-1] + ' ' + self.variableDictionary[variables[len(variables)-1]]
        
        # # Write the cartesian output
        # if len(self.Cartesians) > 0:
            # for i in range(len(self.Cartesians)):
                # for j in range(len(self.Cartesians[i])-1):
                    # cartesianOutput += self.Cartesians[i][j] + ' '
                # cartesianOutput += self.Cartesians[i][len(self.Cartesians[i])-1] + '\n'
        # cartesianOutput = cartesianOutput[:-1]
        
        # root = os.getcwd()
        # if os.path.exists(root + '/zmatFiles'):
            # shutil.rmtree(root + '/zmatFiles')
        
        # os.mkdir('zmatFiles')
        # os.chdir('zmatFiles')
        
        
        # with open('atomList','w') as file:
            # file.write(atomListOutput)
        # with open('bondIndices','w') as file:
            # file.write(bondIndicesOutput)
        # with open('bondVariables','w') as file:
            # file.write(bondVariablesOutput)
        # with open('angleIndices','w') as file:
            # file.write(angleIndicesOutput)
        # with open('angleVariables','w') as file:
            # file.write(angleVariablesOutput)
        # with open('torsionIndices','w') as file:
            # file.write(torsionIndicesOutput)
        # with open('torsionVariables','w') as file:
            # file.write(torsionVariablesOutput)
        # with open('variableDictionary','w') as file:
            # file.write(variableDictionaryOutput)
        # with open('Cartesians','w') as file:
            # file.write(cartesianOutput)






