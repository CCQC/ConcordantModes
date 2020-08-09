import numpy as np
import os
import shutil
import re

root = os.getcwd()
packagepath = os.path.realpath(__file__)
packagepath = packagepath[:-len('/GenLoad.py')]

# Define some regexes
intcRegex = re.compile(r'INTC\[\]')
dotRegex = re.compile(r'Dot')

os.chdir(packagepath + '/../templates')

with open('LoadTemplate.wls','r') as file:
    data = file.readlines()

os.chdir(root)

os.chdir('../zmatFiles')

# Read in ZMAT data
with open('atomList','r') as file:
    atoms = file.read()
with open('bondIndices','r') as file:
    bondIndices = file.read()
with open('bondVariables','r') as file:
    bondVariables = file.read()
with open('angleIndices','r') as file:
    angleIndices = file.read()
with open('angleVariables','r') as file:
    angleVariables = file.read()
with open('torsionIndices','r') as file:
    torsionIndices = file.read()
with open('torsionVariables','r') as file:
    torsionVariables = file.read()
with open('variableDictionary','r') as file:
    variableDictionary = file.read()
with open('Cartesians','r') as file:
    Cartesians = file.read()

os.chdir('..')

with open('disp_option','r') as file:
    disp_file = file.read()

os.chdir(root)


# Split the data into arrays
varDict = {}
disp_data     = disp_file.split('\n')
atoms         = atoms.split('\n')
bondVariables = bondVariables.split('\n')
bondIndices   = bondIndices.split('\n')
for i in range(len(bondIndices)):
    bondIndices[i] = bondIndices[i].split(' ')
angleVariables = angleVariables.split('\n')
angleIndices   = angleIndices.split('\n')
for i in range(len(angleIndices)):
    angleIndices[i] = angleIndices[i].split(' ')
torsionVariables = torsionVariables.split('\n')
torsionIndices   = torsionIndices.split('\n')
for i in range(len(torsionIndices)):
    torsionIndices[i] = torsionIndices[i].split(' ')
variableDictionary = variableDictionary.split('\n')
for i in range(len(variableDictionary)):
    variableDictionary[i] = variableDictionary[i].split(' ')
for i in range(len(variableDictionary)):
    varDict[variableDictionary[i][0]] = float(variableDictionary[i][1])

variables = bondVariables.copy()

if angleIndices[0][0] != '':
    for i in angleVariables:
        variables.append(i)
if torsionIndices[0][0] != '':
    for i in torsionVariables:
        variables.append(i)

Cartesians = Cartesians.split('\n')
for i in range(len(Cartesians)):
    Cartesians[i] = Cartesians[i].split(' ')

# print(bondIndices)
# print(angleIndices)
# print(torsionIndices)

# Fine, I'll hard code the masses too...
massesOutputString = 'mass = {'
for i in atoms:
    massesOutputString += 'm' + i + ','
massesOutputString = massesOutputString[:-1]
massesOutputString += '};\n'
data[12] = massesOutputString

print("disp_data: ")
print(disp_data)
# Option for the disp length
dispOutputString = 'rdisp = '
dispOutputString += str(disp_data[0]) + ';\n'
data[26] = dispOutputString
dispOutputString = 'adisp = '
dispOutputString += str(disp_data[1]) + ';\n'
data[27] = dispOutputString

dispsOutputString = 'rr = {'

for i in bondVariables:
    dispsOutputString += 'rdisp,'
if angleIndices[0][0] != '':
    for i in angleVariables:
        dispsOutputString += 'adisp,'
if torsionIndices[0][0] != '':
    for i in torsionVariables:
        dispsOutputString += 'adisp,'
dispsOutputString = dispsOutputString[:-1]
dispsOutputString += '};'
data[29] = dispsOutputString

clearArray = []
clearString = ''
clearString += 'Clear[INTC'

for i in range(len(atoms)):
    clearString += ', x' + str(i+1)
    clearString += ', y' + str(i+1)
    clearString += ', z' + str(i+1)
    if i%6 == 0 and i != 0:
        clearString += '\n'
        clearArray.append(clearString)
        clearString = '  '
clearString += '];\n'
clearArray.append(clearString)
clearArray.reverse()

for i in clearArray:
    data.insert(35,i)
for i in range(len(data)):
    if re.search(intcRegex,data[i]):
        index = i

intcArray = []
intcString = 'INTC['
for i in range(len(atoms)-1):
    intcString += 'x' + str(i+1) + '_, '
    intcString += 'y' + str(i+1) + '_, '
    intcString += 'z' + str(i+1) + '_, '
    if i%6 == 0 and i != 0:
        intcString += '\n'
        intcArray.append(intcString)
        intcString = '   '
intcString += 'x' + str(len(atoms)) + '_, '
intcString += 'y' + str(len(atoms)) + '_, '
intcString += 'z' + str(len(atoms)) + '_] =\n'
intcArray.append(intcString)
intcArray.reverse()

data[index] = '\n'
for i in intcArray:
    data.insert(index+1,i)
for i in range(len(data)):
    if re.search(dotRegex,data[i]):
        index = i

dotArray = []
dotString = '   {'
# Radii from Cartesians equations
for i in range(len(bondIndices)-1):
    dotString += 'Rab[{x' + str(bondIndices[i][0])
    dotString += ',y' + str(bondIndices[i][0])
    dotString += ',z' + str(bondIndices[i][0])
    dotString += '},{'
    dotString += 'x' + str(bondIndices[i][1])
    dotString += ',y' + str(bondIndices[i][1])
    dotString += ',z' + str(bondIndices[i][1])
    dotString += '}],\n'
    dotArray.append(dotString)
    dotString = '    '
dotString += 'Rab[{x' + str(bondIndices[len(bondIndices)-1][0])
dotString += ',y' + str(bondIndices[len(bondIndices)-1][0])
dotString += ',z' + str(bondIndices[len(bondIndices)-1][0])
dotString += '},{'
dotString += 'x' + str(bondIndices[len(bondIndices)-1][1])
dotString += ',y' + str(bondIndices[len(bondIndices)-1][1])
dotString += ',z' + str(bondIndices[len(bondIndices)-1][1])
dotString += '}]'
if angleIndices[0][0] == '':
    dotString += '\n'
    dotArray.append(dotString)
else:
    dotString += ',\n'
    dotArray.append(dotString)
    dotString = '    '

# Angles from Cartesians equations
if angleIndices[0][0] != '':
    for i in range(len(angleIndices)-1):
        dotString += 'ArcCos[Cos\[Theta]abc[{x' + str(angleIndices[i][0])
        dotString += ',y' + str(angleIndices[i][0])
        dotString += ',z' + str(angleIndices[i][0])
        dotString += '},{x' + str(angleIndices[i][1])
        dotString += ',y' + str(angleIndices[i][1])
        dotString += ',z' + str(angleIndices[i][1])
        dotString += '},{x' + str(angleIndices[i][2])
        dotString += ',y' + str(angleIndices[i][2])
        dotString += ',z' + str(angleIndices[i][2])
        dotString += '}]],\n'
        dotArray.append(dotString)
        dotString = '    '
    dotString += 'ArcCos[Cos\[Theta]abc[{x' + str(angleIndices[len(angleIndices)-1][0])
    dotString += ',y' + str(angleIndices[len(angleIndices)-1][0])
    dotString += ',z' + str(angleIndices[len(angleIndices)-1][0])
    dotString += '},{x' + str(angleIndices[len(angleIndices)-1][1])
    dotString += ',y' + str(angleIndices[len(angleIndices)-1][1])
    dotString += ',z' + str(angleIndices[len(angleIndices)-1][1])
    dotString += '},{x' + str(angleIndices[len(angleIndices)-1][2])
    dotString += ',y' + str(angleIndices[len(angleIndices)-1][2])
    dotString += ',z' + str(angleIndices[len(angleIndices)-1][2])
    dotString += '}]]'
    if torsionIndices[0][0] == '':
        dotString += '\n'
        dotArray.append(dotString)
    else:
        dotString += ',\n'
        dotArray.append(dotString)
        dotString = '    '

# Torsion angles from Cartesians equations
if torsionIndices[0][0] != '':
    for i in range(len(torsionIndices)-1):
        dotString += 'ArcTan[Tan\[Tau]abcd[{x' + str(torsionIndices[i][0])
        dotString += ',y' + str(torsionIndices[i][0])
        dotString += ',z' + str(torsionIndices[i][0])
        dotString += '},{x' + str(torsionIndices[i][1])
        dotString += ',y' + str(torsionIndices[i][1])
        dotString += ',z' + str(torsionIndices[i][1])
        dotString += '},{x' + str(torsionIndices[i][2])
        dotString += ',y' + str(torsionIndices[i][2])
        dotString += ',z' + str(torsionIndices[i][2])
        dotString += '},{x' + str(torsionIndices[i][3])
        dotString += ',y' + str(torsionIndices[i][3])
        dotString += ',z' + str(torsionIndices[i][3])
        dotString += '}]]'
        Condition1 = varDict[torsionVariables[i]] > np.pi/2 and varDict[torsionVariables[i]] <= (3*np.pi)/2
        Condition2 = varDict[torsionVariables[i]] < -np.pi/2 and varDict[torsionVariables[i]] >= -(3*np.pi)/2
        if Condition1 or Condition2:
            dotString += '+\[Pi],\n'
        else:
            dotString += ',\n'
        dotArray.append(dotString)
        dotString = '    '
    dotString += 'ArcTan[Tan\[Tau]abcd[{x' + str(torsionIndices[len(torsionIndices)-1][0])
    dotString += ',y' + str(torsionIndices[len(torsionIndices)-1][0])
    dotString += ',z' + str(torsionIndices[len(torsionIndices)-1][0])
    dotString += '},{x' + str(torsionIndices[len(torsionIndices)-1][1])
    dotString += ',y' + str(torsionIndices[len(torsionIndices)-1][1])
    dotString += ',z' + str(torsionIndices[len(torsionIndices)-1][1])
    dotString += '},{x' + str(torsionIndices[len(torsionIndices)-1][2])
    dotString += ',y' + str(torsionIndices[len(torsionIndices)-1][2])
    dotString += ',z' + str(torsionIndices[len(torsionIndices)-1][2])
    dotString += '},{x' + str(torsionIndices[len(torsionIndices)-1][3])
    dotString += ',y' + str(torsionIndices[len(torsionIndices)-1][3])
    dotString += ',z' + str(torsionIndices[len(torsionIndices)-1][3])
    dotString += '}]]'
    Condition1 = varDict[torsionVariables[len(torsionIndices)-1]] > np.pi/2 and varDict[torsionVariables[len(torsionIndices)-1]] <= (3*np.pi)/2
    Condition2 = varDict[torsionVariables[len(torsionIndices)-1]] < -np.pi/2 and varDict[torsionVariables[len(torsionIndices)-1]] >= -(3*np.pi)/2
    if Condition1 or Condition2:
        dotString += '+\[Pi]\n'
    else:
        dotString += '\n'
    dotArray.append(dotString)
dotString = '    }];\n'
dotArray.append(dotString) 

dotArray.reverse()
for i in dotArray:
    data.insert(index+1,i)

if os.path.exists(os.getcwd() + '/Load.wls'):
    os.remove(os.getcwd() + '/Load.wls')

with open('Load.wls','w') as file:
    file.writelines(data)

os.system('chmod +x Load.wls')

# Here all of the imported files will be generated

# First, the masses
# massesOutputString = ''
# for i in atoms:
    # massesOutputString += 'm' + i +  ','
# massesOutputString = massesOutputString[:-1]

# with open('masses.csv','w') as file:
    # file.write(massesOutputString)

# Next the variables
variablesOutputString = ''
for i in variables:
    variablesOutputString += str(varDict[i]) + ','
variablesOutputString = variablesOutputString[:-1]

with open('variables.csv','w') as file:
    file.write(variablesOutputString)

# Next the disp variables
# dispsOutputString = ''
# for i in bondVariables:
    # dispsOutputString += 'rdisp,'
# if angleIndices[0][0] != '':
    # for i in angleVariables:
        # dispsOutputString += 'adisp,'
# if torsionIndices[0][0] != '':
    # for i in torsionVariables:
        # dispsOutputString += 'adisp,'
# dispsOutputString = dispsOutputString[:-1]

# with open('disps.csv','w') as file:
    # file.write(dispsOutputString)

# Finally, the cartesians
cartesiansOutputString = ''
for i in range(len(Cartesians)):
    cartesiansOutputString += Cartesians[i][0] + ',' + Cartesians[i][1] + ',' + Cartesians[i][2] + '\n'
cartesiansOutputString = cartesiansOutputString[:-1]

with open('cartesians.csv','w') as file:
    file.write(cartesiansOutputString)

# for i in data:
    # print(i)
