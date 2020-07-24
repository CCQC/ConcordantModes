# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 11:00:01 2015

@author: marcus_bartlett
"""

import numpy as np
from numpy.linalg import inv
m = np.loadtxt("eigen.csv", delimiter=',')
# print(m)


ainv = inv(m)

"""
  Below is a generalized version of Marcus' code that should be able to work for arbitrary-sized matrices.
"""
outputString = ''
for i in range(len(ainv)):
    line = ''
    k = 0
    line += "{:5d}".format(i+1)
    line += "{:4d}".format(1)
    # line += "{:+11.6f}".format(ainv[i][0])
    line += "{:+14.9f}".format(ainv[i][0])
    for j in range(len(ainv)-1):
        width = (k+1)*76
        if len(line)<(width-1):
            line += "{:4d}".format(j+2)
            # line += "{:+11.6f}".format(ainv[i][j+1])
            line += "{:+14.9f}".format(ainv[i][j+1])
        else:
            line += "\n"
            k += 1
            line += "{:5d}".format(i+1)
            line += "{:4d}".format(j+2)
            # line += "{:+11.6f}".format(ainv[i][j+1])
            line += "{:+14.9f}".format(ainv[i][j+1])
    outputString += line + '\n'
    # print(line)

with open('SymAdaptInts.dat','w') as file:
    file.write(outputString)

"""
 Below is old code written by Marcus. This is the template for properly formatting a 30x30 symmetry
 internal coordinate matrix for INTDER from a 30 degrees of freedom system from its internal coordinate vibrational eigenvalue matrix.
"""

# for i in range(len(ainv)):
    # print("{:5d}".format(i+1) + "{:4d}".format(1) + "{:+11.6f}".format(ainv[i][0]) + "{:7d}".format(2) + "{:+11.6f}".format(ainv[i][1])+ "{:7d}".format(3) + "{:+11.6f}".format(ainv[i][3]))
    # print("{:5d}".format(i+1) + "{:4d}".format(5) + "{:+11.6f}".format(ainv[i][4]) + "{:7d}".format(6) + "{:+11.6f}".format(ainv[i][5])+ "{:7d}".format(7) + "{:+11.6f}".format(ainv[i][7]))
    # print("{:5d}".format(i+1) + "{:4d}".format(9) + "{:+11.6f}".format(ainv[i][8]) + "{:7d}".format(10) + "{:+11.6f}".format(ainv[i][9])+ "{:7d}".format(11) + "{:+11.1.6f}".format(ainv[i][11]))
    # print("{:5d}".format(i+1) + "{:4d}".format(13) + "{:+11.6f}".format(ainv[i][12]) + "{:7d}".format(14) + "{:+11.6f}".format(ainv[i][13])+ "{:7d}".format(15) + "{:+:+11.6f}".format(ainv[i][15]))
    # print("{:5d}".format(i+1) + "{:4d}".format(17) + "{:+11.6f}".format(ainv[i][16]) + "{:7d}".format(18) + "{:+11.6f}".format(ainv[i][17])+ "{:7d}".format(19) + "{:+:+11.6f}".format(ainv[i][19]))
    # print("{:5d}".format(i+1) + "{:4d}".format(21) + "{:+11.6f}".format(ainv[i][20]) + "{:7d}".format(22) + "{:+11.6f}".format(ainv[i][21])+ "{:7d}".format(23) + "{:+:+11.6f}".format(ainv[i][23]))
    # print("{:5d}".format(i+1) + "{:4d}".format(25) + "{:+11.6f}".format(ainv[i][24]) + "{:7d}".format(26) + "{:+11.6f}".format(ainv[i][25])+ "{:7d}".format(27) + "{:+:+11.6f}".format(ainv[i][27]))
    # print("{:5d}".format(i+1) + "{:4d}".format(29) + "{:+11.6f}".format(ainv[i][28]) + "{:7d}".format(30) + "{:+11.6f}".format(ainv[i][29]))
