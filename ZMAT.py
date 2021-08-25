import numpy as np
import os
import shutil
import re
from . import masses
from ConcordantModes.int2cart import int2cart
from ConcordantModes.TransDisp import TransDisp

class ZMAT(object):
    def __init__(self,options):
        self.amu_elMass = 5.48579909065*(10**(-4))
        self.dispTol    = 1.0e-14
        self.options    = options
        self.Bohr_Ang   = 0.529177210903

    def run(self):
        # Define some regexes
        zmatBeginRegex  = re.compile(r"ZMAT begin")
        zmatEndRegex    = re.compile(r"ZMAT end")
       
        """ ZMAT regexes"""
        firstAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s*\n")
        secondAtomRegex = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s*\n")
        thirdAtomRegex  = re.compile("^\s*([A-Z][a-z]*)\s+(\d+)\s+(\d+)\s*\n")
        fullAtomRegex   = re.compile(
                "^\s*([A-Z][a-z]*)\s+(\d+)\s+(\d+)\s+(\d+)\s*\n")
        """ Custom int coord regexes """
        bondRegex    = re.compile("^\s*(\d+)\s+(\d+)\s*\n")
        angleRegex   = re.compile("^\s*(\d+)\s+(\d+)\s+(\d+)\s*\n")
        torsionRegex = re.compile("^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*T\s*\n")
        oopRegex     = re.compile("^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*O\s*\n")
        linRegex     = re.compile("^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*L\s*\n")
        
        """ Cartesian regexes"""
        cartBeginRegex     = re.compile(r"cart begin")
        cartBeginRegex     = re.compile(r"cart begin")
        cartEndRegex       = re.compile(r"cart end")
        s = "[A-Z][A-Za-z]*\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\n"
        cartesianRegex     = re.compile(s)
        s = "([A-Z][A-Za-z]*)\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n"
        cartesianAtomRegex = re.compile(s)
        dividerRegex       = re.compile("^\s*\-\-\-\s*\n")
        
        """ Read in the ZMAT file """
        with open("zmat",'r') as file:
            output = file.readlines()

        """
            Read in the input cartesian coordinates
        """
        self.CartesiansInit  = []
        self.CartesiansFinal = []
        cartRange            = []
        self.atomList        = []
        
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
        
        for i in range(len(cartOutputInit)):
            if re.search(cartesianRegex,cartOutputInit[i]):
                temp = re.findall(cartesianRegex,cartOutputInit[i])
                atom = re.findall(cartesianAtomRegex,cartOutputInit[i])
                self.CartesiansInit.append(temp[0])
                self.atomList.append(atom[0])
        self.CartesiansInit = np.array(self.CartesiansInit).astype(float)

        if divider:
            for i in range(len(cartOutputFinal)):
                if re.search(cartesianRegex,cartOutputFinal[i]):
                    temp = re.findall(cartesianRegex,cartOutputFinal[i])
                    self.CartesiansFinal.append(temp[0])
            self.CartesiansFinal = np.array(self.CartesiansFinal).astype(float)
        else:
            self.CartesiansFinal = self.CartesiansInit.copy()
       
        if self.options.cartCoords.upper() == "ANGSTROM":
            self.CartesiansInit /= self.Bohr_Ang
            self.CartesiansFinal /= self.Bohr_Ang

        """ Slice out the ZMAT from the input """
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

        """ Initialize necessary lists """
        self.bondIndices             = []
        self.bondVariables           = []
        self.angleIndices            = []
        self.angleVariables          = []
        self.torsionIndices          = []
        self.torsionVariables        = []
        self.oopIndices              = []
        self.oopVariables            = []
        self.linIndices              = []
        self.linVariables            = []
        self.variableDictionaryInit  = {}
        self.variableDictionaryFinal = {}

        count = 0
        if self.options.coords.upper() == "ZMAT":
            """ This code reaps ZMAT data"""
            for i in range(len(zmatOutput)):
                """ This case if we are at the first atom of the ZMAT """
                if re.search(firstAtomRegex,zmatOutput[i]) and count < 1:
                    firstIndex = i
                    count += 1
                """ Second atom of the ZMAT, will have one bond term """
                if re.search(secondAtomRegex,zmatOutput[i]):
                    List = re.findall(secondAtomRegex,zmatOutput[i])[0]
                    self.bondIndices.append([str(i-firstIndex+1),List[1]])
                    self.bondVariables.append("R"+str(i-firstIndex))
                """ Third atom of the ZMAT, will have bond and angle term """
                if re.search(thirdAtomRegex,zmatOutput[i]):
                    List = re.findall(thirdAtomRegex,zmatOutput[i])[0]
                    self.bondIndices.append([str(i-firstIndex+1),List[1]])
                    self.bondVariables.append("R"+str(i-firstIndex))
                    self.angleIndices.append(
                            [str(i-firstIndex+1),List[1],List[2]])
                    self.angleVariables.append("A"+str(i-firstIndex))
                """ 
                    All remaining ZMAT atoms, will have bond, angle, and torsion 
                    term 
                """
                if re.search(fullAtomRegex,zmatOutput[i]):
                    List = re.findall(fullAtomRegex,zmatOutput[i])[0]
                    self.bondIndices.append([str(i-firstIndex+1),List[1]])
                    self.bondVariables.append("R"+str(i-firstIndex))
                    self.angleIndices.append(
                            [str(i-firstIndex+1),List[1],List[2]])
                    self.angleVariables.append("A"+str(i-firstIndex))
                    self.torsionIndices.append(
                            [str(i-firstIndex+1),List[1],List[2],List[3]])
                    self.torsionVariables.append("D"+str(i-firstIndex))
        elif self.options.coords.upper() == "REDUNDANT":
            count = 0
            for i in range(len(zmatOutput)):
                if re.search(bondRegex,zmatOutput[i]):
                    count += 1
                    List = re.findall(bondRegex,zmatOutput[i])[0]
                    self.bondIndices.append(List)
                    self.bondVariables.append("R"+str(count))
            self.bondIndices = np.array(self.bondIndices)
            """ Form all possible angles from bonds """
            self.angleIndices = np.array([])
            count = 0
            for i in range(len(self.bondIndices)):
                for j in range(len(self.bondIndices)-i-1):
                    a = np.setdiff1d(
                            self.bondIndices[i],self.bondIndices[i+j+1])
                    b = np.intersect1d(
                            self.bondIndices[i],self.bondIndices[i+j+1])
                    c = np.setdiff1d(
                            self.bondIndices[i+j+1],self.bondIndices[i])
                    if len(a) and len(b) and len(c):
                        d = np.array([a[0],b[0],c[0]])
                        self.angleIndices = np.append(self.angleIndices,d)
                        count += 1
                        self.angleVariables.append("A"+str(count))
            self.angleIndices = self.angleIndices.reshape((-1,3))
            
            """ Form all possible torsions from angles """
            self.torsionIndices = np.array([])
            count = 0
            for i in range(len(self.angleIndices)):
                for j in range(len(self.angleIndices)-i-1):
                    a = np.setdiff1d(
                            self.angleIndices[i],self.angleIndices[i+j+1])
                    b = np.intersect1d(
                            self.angleIndices[i],self.angleIndices[i+j+1])
                    c = np.setdiff1d(
                            self.angleIndices[i+j+1],self.angleIndices[i])
                    if len(a) and len(b)==2 and len(c):
                        d = np.array([a[0],b[0],b[1],c[0]])
                        self.torsionIndices = np.append(self.torsionIndices,d)
                        count += 1
                        self.torsionVariables.append("D"+str(count))
            self.torsionIndices = self.torsionIndices.reshape((-1,4))

        elif self.options.coords.upper() == "CUSTOM":
            """ 
                This option will allow the user to specify a custom array of 
                internal coordinates. 
            """
            Sum = 0
            blank = 0
            for i in range(len(zmatOutput)):
                if re.search(bondRegex,zmatOutput[i]):
                    List = re.findall(bondRegex,zmatOutput[i])[0]
                    self.bondIndices.append(List)
                    self.bondVariables.append(
                            "R"+str(i+1-Sum+len(self.bondVariables)))
                elif re.search(angleRegex,zmatOutput[i]):
                    List = re.findall(angleRegex,zmatOutput[i])[0]
                    self.angleIndices.append(List)
                    self.angleVariables.append(
                            "A"+str(i+1-Sum+len(self.angleVariables)))
                elif re.search(torsionRegex,zmatOutput[i]):
                    List = re.findall(torsionRegex,zmatOutput[i])[0]
                    self.torsionIndices.append(List)
                    self.torsionVariables.append(
                            "D"+str(i+1-Sum+len(self.torsionVariables)))
                elif re.search(oopRegex,zmatOutput[i]):
                    List = re.findall(oopRegex,zmatOutput[i])[0]
                    self.oopIndices.append(List)
                    self.oopVariables.append(
                            "O"+str(i+1-Sum+len(self.oopVariables)))
                elif re.search(linRegex,zmatOutput[i]):
                    List = re.findall(linRegex,zmatOutput[i])[0]
                    self.linIndices.append(List)
                    self.linVariables.append(
                            "L"+str(i+1-Sum+len(self.linVariables)))
                else:
                    blank += 1
                Sum = len(self.bondVariables) + len(self.angleVariables) \
                    + len(self.torsionVariables) + len(self.oopVariables) \
                    + len(self.linVariables) + blank
        
        """
            This code utilizes the INTC function from the TransDisp module to 
            calculate the initial variable values from the cartesian 
            coordinates.
        """
        indices = []
        transdisp = TransDisp(1,self,1,1,False,self.dispTol,np.array([]),self.options,indices)
        I = np.eye(len(self.bondIndices)+len(self.angleIndices)+len(self.torsionIndices)+len(self.oopIndices))
        variables1 = transdisp.INTC(self.CartesiansInit,I,I)
        variables2 = transdisp.INTC(self.CartesiansFinal,I,I)
        for i in range(
                len(self.angleIndices) + len(self.torsionIndices) \
                + len(self.oopIndices) + len(self.linIndices)):
            variables1[len(self.bondIndices)+i] *= 180./np.pi
            variables2[len(self.bondIndices)+i] *= 180./np.pi
        Variables = np.append(
                self.bondVariables,np.append(self.angleVariables,
                    self.torsionVariables))
        if len(self.oopVariables):
            Variables = np.append(Variables,self.oopVariables)
        if len(self.linVariables):
            Variables = np.append(Variables,self.linVariables)

        """
            This code is useful for checking the Redundant coordinate 
            generation process.
        """
        if self.options.coords.upper() == "REDUNDANT":
            indices = []
            for i in range(len(self.bondIndices)):
                indices.append(self.bondIndices[i].tolist())
            for i in range(len(self.angleIndices)):
                indices.append(self.angleIndices[i].tolist())
            for i in range(len(self.torsionIndices)):
                indices.append(self.torsionIndices[i].tolist())
            # print(indices)

        for i in range(len(variables1)):
            self.variableDictionaryInit[Variables[i]] = variables1[i]

        if divider:
            for i in range(len(variables2)):
                self.variableDictionaryFinal[Variables[i]] = variables2[i]
        else:
            self.variableDictionaryFinal = self.variableDictionaryInit.copy()
        """ 
            And now we must temper the torsion angles! For consistency's sake 
            we will force them to lie between -90 deg and +270 deg.

            *** Out of plane and linear bends will have to be added
        """
        
        """ Handle Variable lists separately. First the INIT: """
        for i in range(len(self.torsionVariables)):
            Condition1 = float(
                    self.variableDictionaryInit[
                        self.torsionVariables[i]]) <= -90.
            Condition2 = float(
                    self.variableDictionaryInit[
                        self.torsionVariables[i]]) >= 270.
            buff = np.floor(abs(float(self.variableDictionaryInit[
                self.torsionVariables[i]]))/360)
            if Condition1:
                self.variableDictionaryInit[self.torsionVariables[i]] = \
                    float(self.variableDictionaryInit[self.torsionVariables[i]])
                self.variableDictionaryInit[self.torsionVariables[i]] += \
                    360.*buff
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
            The masses are assigned to the respective atom from the masses.py file
        """
        self.masses = [masses.get_mass(label) for label in self.atomList]
        for i in range(len(self.masses)):
            self.masses[i] = self.masses[i]/self.amu_elMass
       
        """
            Print off the internal coordinate and its value in Bohr/Degree
        """
        print("Initial Geometric Internal Coordinate Values:")
        for i in range(len(Variables)):
            print(Variables[i] + " = " + str(self.variableDictionaryInit[Variables[i]]))
            # if self.options.coords.upper() == "REDUNDANT":
                # print(indices[i])
        print("Final Geometric Internal Coordinate Values:")
        for i in range(len(Variables)):
            print(Variables[i] + " = " + str(self.variableDictionaryFinal[Variables[i]]))
            # if self.options.coords.upper() == "REDUNDANT":
                # print(indices[i])
        print("Final - Initial Geometric Internal Coordinate Values:")
        for i in range(len(Variables)):
            print(Variables[i] + " = " + str(self.variableDictionaryFinal[Variables[i]] - self.variableDictionaryInit[Variables[i]]))
            # if self.options.coords.upper() == "REDUNDANT":
                # print(indices[i])
        """
            Calculate Cartesians using ZMATs: Sadly this will have to go on the backburner.
            The cartesians must match those used to generate the cartesian force constants or you're gonna have a bad time.
        """
        # cartInit = int2cart(self,self.variableDictionaryInit)
        # cartFinal = int2cart(self,self.variableDictionaryFinal)
        # cartInit.run()
        # cartFinal.run()
        # self.CartesiansInit = cartInit.Carts
        # self.CartesiansFinal = cartFinal.Carts
        # print(self.CartesiansInit)
        # print(self.CartesiansFinal)
        
        
        # print(self.CartesiansInit)
        # print(self.CartesiansFinal)
        # print(self.bondVariables)
        # print(self.torsionIndices)
