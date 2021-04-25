import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power

"""
    This class will be used to print out total energy distributions
    from internal coordinate eigenvalue matrices. If redundant
    coordinates are utilized, projection matrices will be stored here.
"""

class TED(object):
    def __init__(self, Proj, Proj_inv, zmat):
        self.Proj     = Proj
        self.Proj_inv = Proj_inv
        self.zmat     = zmat

    def run(self, eigs, Freq):
        """
            Need to properly form this TED equation.
        """
        # TED = np.multiply(np.dot(self.L_proj,eigs),np.dot(self.R_proj,np.transpose(LA.inv(eigs))))*100
        self.table_print(Freq,TED)

    def table_print(self, Freq, TED):
        tableOutput = "{:>26s}".format("Frequency #: ")
        for i in range(len(Freq)):
            tableOutput += "{:8d}".format(i+1)
        tableOutput += '\n'
        tableOutput += "--------------------------"
        for i in range(len(Freq)):
            tableOutput += "--------"
        tableOutput += '\n'
        tableOutput += "{:>26s}".format("Frequency: ")
        for i in range(len(Freq)):
            tableOutput += " " + "{:7.1f}".format(Freq[i])
        tableOutput += '\n'
        tableOutput += "--------------------------"
        for i in range(len(Freq)):
            tableOutput += "--------"
        tableOutput += '\n'
        """ Modify table rows so that internal coord labels have a fixed length. """
        for i in range(len(TED)):
            if i < len(self.zmat.bondIndices):
                tableOutput += "{:10s}".format(" ") + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.bondIndices[i][0])-1]) + str(self.zmat.bondIndices[i][0])) + " " \
                        + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.bondIndices[i][1])-1]) + str(self.zmat.bondIndices[i][1])) + " STRE: " 
            elif i < len(self.zmat.bondIndices) + len(self.zmat.angleIndices):
                k = i - len(self.zmat.bondIndices)
                tableOutput += "{:5s}".format(" ") + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.angleIndices[k][0])-1]) + str(self.zmat.angleIndices[k][0])) + " " \
                        + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.angleIndices[k][1])-1]) + str(self.zmat.angleIndices[k][1])) + " " \
                        + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.angleIndices[k][2])-1]) + str(self.zmat.angleIndices[k][2])) + " BEND: " 
            elif i < len(self.zmat.bondIndices) + len(self.zmat.angleIndices) + len(self.zmat.torsionIndices):
                k = i - len(self.zmat.bondIndices) - len(self.zmat.angleIndices)
                tableOutput += "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][0])-1]) + str(self.zmat.torsionIndices[k][0])) + " " \
                        + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][1])-1]) + str(self.zmat.torsionIndices[k][1])) + " " \
                        + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][2])-1]) + str(self.zmat.torsionIndices[k][2])) + " " \
                        + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][3])-1]) + str(self.zmat.torsionIndices[k][3])) + " TORS: " 
            for j in range(len(TED[0])):
                tableOutput += "{:8.1f}".format(TED[i][j])
            tableOutput += '\n'
