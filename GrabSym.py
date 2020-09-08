import re
import os
import numpy as np

class GrabSym(object):
    def __init__(self):
        # Define regexes
        self.symBegin = re.compile("\s*VALUES OF SYMMETRY INTERNAL COORDINATES \(ANG. OR RAD.\)")
        self.symEnd   = re.compile(r"DETERMINANT OF")
        self.symFind  = re.compile("\s*\d+\s*(-?\d+\.\d+)")

    def run(self):
        # First, read in file
        with open('intder.out','r') as file:
            data = file.readlines()
        
        # These lines pull out the proper sym coordinates to search within the intder.out file
        indices = []
        for i in range(len(data)):
            if re.search(self.symBegin,data[i]):
                indices.append(i)
        endIndex = []
        for i in range(len(data)):
            if re.search(self.symEnd,data[i]):
                endIndex.append(i)
        
        index = indices[0]
        lines = []
        for i in range(endIndex[0] - index):
            lines.append(data[i+index])
        
        # Search for the sym coordinates
        symCoords = np.array([])
        for i in lines:
            symCoords = np.append(symCoords,re.findall(self.symFind,i))
        
        # Construct outputString in .csv format
        outputString = ''
        for i in symCoords:
            outputString += i + ','
        outputString = outputString[:-1]
        
        with open('symVariables.csv','w') as file:
            file.write(outputString)
