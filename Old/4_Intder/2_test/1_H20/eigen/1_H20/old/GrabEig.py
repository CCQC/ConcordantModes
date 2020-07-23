import re
import os
import numpy as np

# Define complicated regex
EigRegex = re.compile("VIBRATIONAL\sFREQUENCIES\s\(CM\-1\)\sAND\sEIGENVECTORS\n\n\n\s----------------------------------------\n(?:\s*\d+)+\s*\n\sCOORDINATE(?:\s*-?\d+\.\d+)+\n\s----------------------------------------\n(?:\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*\n)+") 
EigReg   = re.compile("VIBRATIONAL\sFREQUENCIES\s\(CM\-1\)\sAND\sEIGENVECTORS") 
NARegex  = re.compile("^\s+(\d)")
MoreEigs = re.compile("^\s+\d+\s+(-?\d\.\d+)\s+(-?\d\.\d+)\s+(-?\d\.\d+)")
# EigRegexDotAll = re.compile("VIBRATIONAL\sFREQUENCIES\s\(CM\-1\)\sAND\sEIGENVECTORS\n\n\n\s----------------------------------------\n(?:\s*\d+)+\s*\n\sCOORDINATE(?:\s*-?\d+\.\d+)+\n\s----------------------------------------\n(.+) ----------------------------------------",re.DOTALL) 
# EigRegexDotAll = re.compile("VIBRATIONAL\sFREQUENCIES\s\(CM\-1\)\sAND\sEIGENVECTORS\n\n\n\s----------------------------------------\n(?:\s*\d+)+\s*\n\sCOORDINATE(?:\s*-?\d+\.\d+)+\n\s----------------------------------------\n(.+)\s+(?:\-)+\n",re.DOTALL) 


# First, read in file
with open('intder.out','r')  as file:
    data = file.readlines()

# This code is highly ad hoc to determine the dimensionality of the eigenvector matrix. Should be DegFree x DegFree. However there will need to be a check for linearity.
NA = int(re.findall(NARegex,data[11])[0])
DegFree = 3*NA - 6

indices = []
for i in range(len(data)):
    if re.search(EigReg,data[i]):
        indices.append(i)

index = indices[0]
lines = []
for i in range(DegFree):
    lines.append(data[index+i+7])

eigs = []
for i in lines:
    eigs.append(re.findall(MoreEigs,i))

outputString = ''

for i in eigs:
    for j in i[0]:
        outputString += j
        outputString += ','
    outputString = outputString[:-1]
    outputString += '\n'

with open('eigen.csv','w') as file:
    file.write( outputString )



