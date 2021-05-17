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
    def __init__(self, Proj, zmat):
        self.Proj     = Proj
        self.zmat     = zmat

    def run(self, eigs, Freq):
        """
            Need to properly form this TED equation.
        """
        projEigs = np.dot(self.Proj,eigs)
        projEigsInv = LA.pinv(projEigs)
        TED = np.multiply(projEigs,projEigsInv.T)*100
        self.table_print(Freq,TED)

    """
        I will likely want to add two other tables. One that displays the TED of the modes in terms of symmetry adapted coordinates. And another that displays the TED of internals within each symmetry adapted coordinate.
    """
    def table_print(self, Freq, TED):
        if len(Freq) != len(TED[0]):
            print("Something has gone terribly wrong. Your Total Energy distribution should have the same number of columns as frequencies, but it does not.")
            print("TED columns: " + str(len(TED[0])))
            print("# Frequencies: " + str(len(Freq)))
            raise RuntimeError
        div = 15
        a = len(Freq)//div
        for i in range(a):
            print(self.subTable(Freq[i*div:(i+1)*div], TED, i, div))
        print(self.subTable(Freq[a*div:], TED, a, div))

    def subTable(self, Freq, TED, intDiv, div):
        tableOutput = "{:>26s}".format("Frequency #: ")
        n = len(Freq)
        for l in range(n):
            tableOutput += "{:8d}".format(l + intDiv*div + 1)
        tableOutput += '\n'
        tableOutput += "--------------------------"
        for l in range(n):
            tableOutput += "--------"
        tableOutput += '\n'
        tableOutput += "{:>26s}".format("Frequency: ")
        for l in range(n):
            tableOutput += " " + "{:7.1f}".format(Freq[l])
        tableOutput += '\n'
        tableOutput += "--------------------------"
        for l in range(n):
            tableOutput += "--------"
        tableOutput += '\n'
        for i in range(len(TED)):
            for j in range(len(Freq)):
                if i < len(self.zmat.bondIndices) and j == 0:
                    tableOutput += "{:10s}".format(" ") + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.bondIndices[i][0])-1]) + str(self.zmat.bondIndices[i][0])) + " " \
                            + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.bondIndices[i][1])-1]) + str(self.zmat.bondIndices[i][1])) + " STRE: " 
                elif i < len(self.zmat.bondIndices) + len(self.zmat.angleIndices) and j == 0:
                    k = i - len(self.zmat.bondIndices)
                    tableOutput += "{:5s}".format(" ") + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.angleIndices[k][0])-1]) + str(self.zmat.angleIndices[k][0])) + " " \
                            + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.angleIndices[k][1])-1]) + str(self.zmat.angleIndices[k][1])) + " " \
                            + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.angleIndices[k][2])-1]) + str(self.zmat.angleIndices[k][2])) + " BEND: " 
                elif i < len(self.zmat.bondIndices) + len(self.zmat.angleIndices) + len(self.zmat.torsionIndices) and j == 0:
                    k = i - len(self.zmat.bondIndices) - len(self.zmat.angleIndices)
                    tableOutput += "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][0])-1]) + str(self.zmat.torsionIndices[k][0])) + " " \
                            + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][1])-1]) + str(self.zmat.torsionIndices[k][1])) + " " \
                            + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][2])-1]) + str(self.zmat.torsionIndices[k][2])) + " " \
                            + "{:4s}".format(str(self.zmat.atomList[int(self.zmat.torsionIndices[k][3])-1]) + str(self.zmat.torsionIndices[k][3])) + " TORS: " 
                tableOutput += "{:8.1f}".format(TED[i][j])
            tableOutput += '\n'
        
        return tableOutput
