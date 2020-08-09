import re
import os
import numpy as np

# Define complicated regex
eigReg     = re.compile("VIBRATIONAL\sFREQUENCIES\s\(CM\-1\)\sAND\sEIGENVECTORS")
endEigReg  = re.compile("TOTAL\sENERGY\sDISTRIBUTIONS")
rawEigsReg = re.compile("^\s*\d+\s*((?:-?\d+\.\d+\s*)+)")
eigsReg    = re.compile("(-?\d+\.\d+)+")
NARegex    = re.compile("^\s+(\d+)")
# EigRegexDotAll = re.compile("VIBRATIONAL\sFREQUENCIES\s\(CM\-1\)\sAND\sEIGENVECTORS\n\n\n\s----------------------------------------\n(?:\s*\d+)+\s*\n\sCOORDINATE(?:\s*-?\d+\.\d+)+\n\s----------------------------------------\n(.+) ----------------------------------------",re.DOTALL) 

# First, read in file
with open('intder.out','r') as file:
    data = file.readlines()

# This code is highly ad hoc to determine the dimensionality of the eigenvector matrix. Should be degFree x degFree. However there will need to be a check for geometric linearity.
NA = int(re.findall(NARegex,data[11])[0])
degFree = 3*NA - 6

# These lines pull out the proper eigenvalue tables to search within the intder.out file
indices = []
for i in range(len(data)):
    if re.search(eigReg,data[i]):
        indices.append(i)
endIndex = []
for i in range(len(data)):
    if re.search(endEigReg,data[i]):
        endIndex.append(i)

index = indices[0]
lines = []
for i in range(endIndex[0] - index):
    lines.append(data[i+index])


# This code pulls out all of the eigenvalues, but in a non-sorted format, and converts the data type from strings to floats
rawEigs = np.array([])
for i in lines:
    rawEigs = np.append(rawEigs,re.findall(rawEigsReg,i))

eigs = np.array([])
for i in rawEigs:
    eigs = np.append(eigs,re.findall(eigsReg,i))

eigs = eigs.astype(float)


# These lines sort the eigenvalues into a proper matrix format. There isn't an easy way to understand what's
# going on here, unfortunately.
sortedEigs = np.zeros((degFree,degFree))
k = 0
l = 0
m = degFree//10
decrement = 10 - degFree%10
increment = 10
if degFree > 9:
    for i in range(degFree):
        for j in range(degFree):
            if j + 1 > m*10:
                l = 1
            sortedEigs[i][j] = eigs[i*(increment - l*decrement) + j + k*(10*(degFree - 1))]
            if (j+1)%10 == 0:
                k += 1
        k = 0
        l = 0
else:
    sortedEigs = np.reshape(eigs,(degFree,degFree))

# These final lines format the eigenvalue matrix into an 'eigen.csv' file that may be further processed by the TED.py script.
outputString = ''

for i in sortedEigs:
    for j in i:
        outputString += str(j)
        outputString += ','
    outputString = outputString[:-1]
    outputString += '\n'

with open('eigen.csv','w') as file:
    file.write( outputString )



