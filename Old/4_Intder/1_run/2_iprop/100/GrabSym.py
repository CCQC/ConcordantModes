import re
import os
import numpy as np

# Define complicated regex
symBegin     = re.compile("\s*VALUES OF SYMMETRY INTERNAL COORDINATES \(ANG. OR RAD.\)")
symEnd     = re.compile(r"DETERMINANT OF")
symFind = re.compile("\s*\d+\s*(-?\d+\.\d+)")

# First, read in file
with open('intder.out','r') as file:
    data = file.readlines()

# These lines pull out the proper sym coordinates to search within the intder.out file
indices = []
for i in range(len(data)):
    if re.search(symBegin,data[i]):
        indices.append(i)
endIndex = []
for i in range(len(data)):
    if re.search(symEnd,data[i]):
        endIndex.append(i)

index = indices[0]
lines = []
for i in range(endIndex[0] - index):
    lines.append(data[i+index])

# Search for the sym coordinates
symCoords = np.array([])
for i in lines:
    symCoords = np.append(symCoords,re.findall(symFind,i))
print(symCoords)

