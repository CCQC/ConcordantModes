import numpy as np
import os
import shutil
import re

class GenLoad(object):
    def __init__(self,rdisp,adisp,zmat):
        self.rdisp = rdisp
        self.adisp = adisp
        self.zmat = zmat

    def run(self):
        root = os.getcwd()
        packagepath = os.path.realpath(__file__)
        packagepath = packagepath[:-len('/GenLoad.py')]
        
        # Define some regexes
        intcRegex = re.compile(r'INTC\[\]')
        dotRegex = re.compile(r'Dot')
        
        os.chdir(packagepath + '/templates')
        
        with open('LoadTemplate.wls','r') as file:
            data = file.readlines()
        
        os.chdir(root)
        
        # Fine, I'll hard code the masses too...
        massesOutputString = 'mass = {'
        for i in self.zmat.atomList:
            massesOutputString += 'm' + i + ','
        massesOutputString = massesOutputString[:-1]
        massesOutputString += '};\n'
        data[12] = massesOutputString
        
        # Option for the disp length
        dispOutputString = 'rdisp = '
        dispOutputString += str(self.rdisp) + ';\n'
        data[26] = dispOutputString
        dispOutputString = 'adisp = '
        dispOutputString += str(self.adisp) + ';\n'
        data[27] = dispOutputString
        
        dispsOutputString = 'rr = {'
        
        for i in self.zmat.bondVariables:
            dispsOutputString += 'rdisp,'
        if len(self.zmat.angleIndices) > 0:
            for i in self.zmat.angleVariables:
                dispsOutputString += 'adisp,'
        if len(self.zmat.torsionIndices) > 0:
            for i in self.zmat.torsionVariables:
                dispsOutputString += 'adisp,'
        dispsOutputString = dispsOutputString[:-1]
        dispsOutputString += '};'
        data[29] = dispsOutputString
        
        clearArray = []
        clearString = ''
        clearString += 'Clear[INTC'
        
        for i in range(len(self.zmat.atomList)):
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
        for i in range(len(self.zmat.atomList)-1):
            intcString += 'x' + str(i+1) + '_, '
            intcString += 'y' + str(i+1) + '_, '
            intcString += 'z' + str(i+1) + '_, '
            if i%6 == 0 and i != 0:
                intcString += '\n'
                intcArray.append(intcString)
                intcString = '   '
        intcString += 'x' + str(len(self.zmat.atomList)) + '_, '
        intcString += 'y' + str(len(self.zmat.atomList)) + '_, '
        intcString += 'z' + str(len(self.zmat.atomList)) + '_] =\n'
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
        for i in range(len(self.zmat.bondIndices)-1):
            dotString += 'Rab[{x' + str(self.zmat.bondIndices[i][0])
            dotString += ',y' + str(self.zmat.bondIndices[i][0])
            dotString += ',z' + str(self.zmat.bondIndices[i][0])
            dotString += '},{'
            dotString += 'x' + str(self.zmat.bondIndices[i][1])
            dotString += ',y' + str(self.zmat.bondIndices[i][1])
            dotString += ',z' + str(self.zmat.bondIndices[i][1])
            dotString += '}],\n'
            dotArray.append(dotString)
            dotString = '    '
        dotString += 'Rab[{x' + str(self.zmat.bondIndices[len(self.zmat.bondIndices)-1][0])
        dotString += ',y' + str(self.zmat.bondIndices[len(self.zmat.bondIndices)-1][0])
        dotString += ',z' + str(self.zmat.bondIndices[len(self.zmat.bondIndices)-1][0])
        dotString += '},{'
        dotString += 'x' + str(self.zmat.bondIndices[len(self.zmat.bondIndices)-1][1])
        dotString += ',y' + str(self.zmat.bondIndices[len(self.zmat.bondIndices)-1][1])
        dotString += ',z' + str(self.zmat.bondIndices[len(self.zmat.bondIndices)-1][1])
        dotString += '}]'
        if len(self.zmat.angleIndices) == 0:
            dotString += '\n'
            dotArray.append(dotString)
        else:
            dotString += ',\n'
            dotArray.append(dotString)
            dotString = '    '
        
        # Angles from Cartesians equations
        if len(self.zmat.angleIndices) > 0:
            for i in range(len(self.zmat.angleIndices)-1):
                dotString += 'ArcCos[Cos\[Theta]abc[{x' + str(self.zmat.angleIndices[i][0])
                dotString += ',y' + str(self.zmat.angleIndices[i][0])
                dotString += ',z' + str(self.zmat.angleIndices[i][0])
                dotString += '},{x' + str(self.zmat.angleIndices[i][1])
                dotString += ',y' + str(self.zmat.angleIndices[i][1])
                dotString += ',z' + str(self.zmat.angleIndices[i][1])
                dotString += '},{x' + str(self.zmat.angleIndices[i][2])
                dotString += ',y' + str(self.zmat.angleIndices[i][2])
                dotString += ',z' + str(self.zmat.angleIndices[i][2])
                dotString += '}]],\n'
                dotArray.append(dotString)
                dotString = '    '
            dotString += 'ArcCos[Cos\[Theta]abc[{x' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][0])
            dotString += ',y' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][0])
            dotString += ',z' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][0])
            dotString += '},{x' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][1])
            dotString += ',y' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][1])
            dotString += ',z' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][1])
            dotString += '},{x' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][2])
            dotString += ',y' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][2])
            dotString += ',z' + str(self.zmat.angleIndices[len(self.zmat.angleIndices)-1][2])
            dotString += '}]]'
            if len(self.zmat.torsionIndices) == 0:
                dotString += '\n'
                dotArray.append(dotString)
            else:
                dotString += ',\n'
                dotArray.append(dotString)
                dotString = '    '
        
        # Torsion angles from Cartesians equations
        if len(self.zmat.torsionIndices) > 0:
            for i in range(len(self.zmat.torsionIndices)-1):
                dotString += 'ArcTan[Tan\[Tau]abcd[{x' + str(self.zmat.torsionIndices[i][0])
                dotString += ',y' + str(self.zmat.torsionIndices[i][0])
                dotString += ',z' + str(self.zmat.torsionIndices[i][0])
                dotString += '},{x' + str(self.zmat.torsionIndices[i][1])
                dotString += ',y' + str(self.zmat.torsionIndices[i][1])
                dotString += ',z' + str(self.zmat.torsionIndices[i][1])
                dotString += '},{x' + str(self.zmat.torsionIndices[i][2])
                dotString += ',y' + str(self.zmat.torsionIndices[i][2])
                dotString += ',z' + str(self.zmat.torsionIndices[i][2])
                dotString += '},{x' + str(self.zmat.torsionIndices[i][3])
                dotString += ',y' + str(self.zmat.torsionIndices[i][3])
                dotString += ',z' + str(self.zmat.torsionIndices[i][3])
                dotString += '}]]'
                Condition1 = float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) > np.pi/2 and float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) <= (3*np.pi)/2
                Condition2 = float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) < -np.pi/2 and float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) >= -(3*np.pi)/2
                if Condition1 or Condition2:
                    dotString += '+\[Pi],\n'
                else:
                    dotString += ',\n'
                dotArray.append(dotString)
                dotString = '    '
            dotString += 'ArcTan[Tan\[Tau]abcd[{x' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][0])
            dotString += ',y' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][0])
            dotString += ',z' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][0])
            dotString += '},{x' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][1])
            dotString += ',y' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][1])
            dotString += ',z' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][1])
            dotString += '},{x' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][2])
            dotString += ',y' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][2])
            dotString += ',z' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][2])
            dotString += '},{x' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][3])
            dotString += ',y' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][3])
            dotString += ',z' + str(self.zmat.torsionIndices[len(self.zmat.torsionIndices)-1][3])
            dotString += '}]]'
            Condition1 = float(self.zmat.variableDictionary[self.zmat.torsionVariables[len(self.zmat.torsionIndices)-1]]) > np.pi/2 and float(self.zmat.variableDictionary[self.zmat.torsionVariables[len(self.zmat.torsionIndices)-1]]) <= (3*np.pi)/2
            Condition2 = float(self.zmat.variableDictionary[self.zmat.torsionVariables[len(self.zmat.torsionIndices)-1]]) < -np.pi/2 and float(self.zmat.variableDictionary[self.zmat.torsionVariables[len(self.zmat.torsionIndices)-1]]) >= -(3*np.pi)/2
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
        
        # Here the imported cartesian files will be generated
        cartesiansOutputString = ''
        for i in range(len(self.zmat.Cartesians)):
            cartesiansOutputString += self.zmat.Cartesians[i][0] + ',' + self.zmat.Cartesians[i][1] + ',' + self.zmat.Cartesians[i][2] + '\n'
        cartesiansOutputString = cartesiansOutputString[:-1]
        
        with open('cartesians.csv','w') as file:
            file.write(cartesiansOutputString)
