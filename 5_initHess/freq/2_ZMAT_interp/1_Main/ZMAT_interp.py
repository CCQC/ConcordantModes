import numpy as np
import os
import shutil
import re


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
atomList           = []
bondIndices        = []
bondVariables      = []
angleIndices       = []
angleVariables     = []
torsionIndices     = []
torsionVariables   = []
variableDictionary = {}
# Now to reap the ZMAT data 
for i in range(len(output)):
    # This case if we are at the first atom of the ZMAT
    if re.search(firstAtomRegex,output[i]):
        atomList.append(re.findall(firstAtomRegex,output[i])[0])
        firstIndex = i
    # Second atom of the ZMAT, will have one bond term
    if re.search(secondAtomRegex,output[i]):
        List = re.findall(secondAtomRegex,output[i])[0]
        atomList.append(List[0])
        bondIndices.append([str(i-firstIndex+1),List[1]])
        bondVariables.append(List[2])
    # Third atom of the ZMAT, will have bond and angle term
    if re.search(thirdAtomRegex,output[i]):
        List = re.findall(thirdAtomRegex,output[i])[0]
        atomList.append(List[0])
        bondIndices.append([str(i-firstIndex+1),List[1]])
        bondVariables.append(List[2])
        angleIndices.append([str(i-firstIndex+1),List[1],List[3]])
        angleVariables.append(List[4])
    # All remaining ZMAT atoms, will have bond, angle, and torsion term
    if re.search(fullAtomRegex,output[i]):
        List = re.findall(fullAtomRegex,output[i])[0]
        atomList.append(List[0])
        bondIndices.append([str(i-firstIndex+1),List[1]])
        bondVariables.append(List[2])
        angleIndices.append([str(i-firstIndex+1),List[1],List[3]])
        angleVariables.append(List[4])
        torsionIndices.append([str(i-firstIndex+1),List[1],List[3],List[5]])
        torsionVariables.append(List[6])
    # Now to reap the variables and their values
    if re.search(variableRegex,output[i]):
        List = re.findall(variableRegex,output[i])[0]
        variableDictionary[List[0]] = List[1]

# Make a cartesian output slice
for i in range(len(tempOutput)):
    if re.search(cartesianRegex,tempOutput[i]):
        cartIndex = i
        break

# Find the cartesian coordinates, in bohr
Cartesians = []
for i in range(len(tempOutput)):
    if re.search(cartesianRegex,tempOutput[i]):
        temp = re.findall(cartesianRegex,tempOutput[i])
        Cartesians.append(temp[0])
    if i > (cartIndex + len(atomList) + 1):
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
if len(atomList) > 0:
    for i in range(len(atomList)-1):
        atomListOutput += atomList[i] + '\n'
    atomListOutput += atomList[len(atomList)-1]

# Write the bond outputs
if len(bondIndices) > 0:
    for i in range(len(bondIndices)-1):
        bondIndicesOutput += bondIndices[i][0] + ' ' + bondIndices[i][1] + '\n'
    bondIndicesOutput += bondIndices[len(bondIndices)-1][0] + ' ' + bondIndices[len(bondIndices)-1][1]
    for i in range(len(bondVariables)-1):
        bondVariablesOutput += bondVariables[i] + '\n'
    bondVariablesOutput += bondVariables[len(bondVariables)-1]

# Write the angle outputs
if len(angleIndices) > 0:
    for i in range(len(angleIndices)-1):
        angleIndicesOutput += angleIndices[i][0] + ' ' + angleIndices[i][1] + ' ' + angleIndices[i][2] + '\n'
    angleIndicesOutput += angleIndices[len(angleIndices)-1][0] + ' ' + angleIndices[len(angleIndices)-1][1] + ' ' + angleIndices[len(angleIndices)-1][2]
    for i in range(len(angleVariables)-1):
        angleVariablesOutput += angleVariables[i] + '\n'
    angleVariablesOutput += angleVariables[len(angleVariables)-1]

# Write the torsion outputs
if len(torsionIndices) > 0:
    for i in range(len(torsionIndices)-1):
        torsionIndicesOutput += torsionIndices[i][0] + ' ' + torsionIndices[i][1] + ' ' + torsionIndices[i][2] + ' ' + torsionIndices[i][3] + '\n'
    torsionIndicesOutput += torsionIndices[len(torsionIndices)-1][0] + ' ' + torsionIndices[len(torsionIndices)-1][1] + ' ' + torsionIndices[len(torsionIndices)-1][2] + ' ' + torsionIndices[len(torsionIndices)-1][3]
    for i in range(len(torsionVariables)-1):
        torsionVariablesOutput += torsionVariables[i] + '\n'
    torsionVariablesOutput += torsionVariables[len(torsionVariables)-1]

# Write out the variables and their values
bondVariables = list(dict.fromkeys(bondVariables))
angleVariables = list(dict.fromkeys(angleVariables))
torsionVariables = list(dict.fromkeys(torsionVariables))
variables = np.array(bondVariables)
variables = np.append(variables,angleVariables)
variables = np.append(variables,torsionVariables)
if len(variableDictionary) > 0:
    if len(variables) > 0:
        for i in range(len(variables)-1):
            variableDictionaryOutput += variables[i] + ' ' + variableDictionary[variables[i]] + '\n'
        variableDictionaryOutput += variables[len(variables)-1] + ' ' + variableDictionary[variables[len(variables)-1]]

# Write the cartesian output
if len(Cartesians) > 0:
    for i in range(len(Cartesians)):
        for j in range(len(Cartesians[i])-1):
            cartesianOutput += Cartesians[i][j] + ' '
        cartesianOutput += Cartesians[i][len(Cartesians[i])-1] + '\n'
cartesianOutput = cartesianOutput[:-1]

with open('atomList','w') as file:
    file.write(atomListOutput)
with open('bondIndices','w') as file:
    file.write(bondIndicesOutput)
with open('bondVariables','w') as file:
    file.write(bondVariablesOutput)
with open('angleIndices','w') as file:
    file.write(angleIndicesOutput)
with open('angleVariables','w') as file:
    file.write(angleVariablesOutput)
with open('torsionIndices','w') as file:
    file.write(torsionIndicesOutput)
with open('torsionVariables','w') as file:
    file.write(torsionVariablesOutput)
with open('variableDictionary','w') as file:
    file.write(variableDictionaryOutput)
with open('Cartesians','w') as file:
    file.write(cartesianOutput)






