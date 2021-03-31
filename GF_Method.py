import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power

"""
    This class will be used as a versatile method to compute GF frequencies!
"""

class GF_Method(object):
    def __init__(self, G, F, tol, projTol, zmat):
        self.G               = G
        self.F               = F
        self.tol             = tol
        self.projTol         = projTol
        self.zmat            = zmat
        self.amu_elMass      = 5.48579909065*(10**(-4))
        self.Hartree_wavenum = 219474.6313708

    def run(self):
        """ Construct the orthogonalizer """
        self.G_O = fractional_matrix_power(self.G,0.5)
        """ Symmetrize F, diagonalize, then backtransform the eigenvectors. """
        self.F_O = np.dot(np.dot(self.G_O,self.F),self.G_O)
        # self.eig_v, self.L_p = LA.eigh(self.F_O)
        self.L_p,self.eig_v,c = LA.svd(self.F_O)
        print("Eigs:")
        print(np.flip(self.eig_v,axis=0))
        # print("EigVecs:")
        # print(self.L_p)
        self.eig_v = np.flip(self.eig_v,axis=0)
        self.L_p = np.flip(self.L_p,axis=1)
        self.L_p[np.abs(self.L_p) < self.tol] = 0
        self.L = np.dot(self.G_O,self.L_p)
        self.L[np.abs(self.L**2) < self.tol] = 0
        self.L = np.real(self.L)
        """ Project out the singular modes """
        self.project()
        print("EigVecs:")
        print(self.L)
        print(np.dot(self.L_p,np.transpose(self.L_p)))
        """ Compute the frequencies by the square root of the eigenvalues. """
        self.Freq = np.sqrt(self.eig_v, dtype=np.complex)
        """ Filter for imaginary modes. """
        for i in range(len(self.Freq)):
            if np.real(self.Freq[i]) > 0.:
                self.Freq[i] = np.real(self.Freq[i])
            else:
                self.Freq[i] = -np.imag(self.Freq[i])
        self.Freq = self.Freq.astype(float)
        """ Convert from Hartrees to wavenumbers. """
        self.Freq *= self.Hartree_wavenum
        # print('Frequencies:')
        # print(self.Freq)
        for i in range(len(self.Freq)):
            print("Frequency #" + "{:3d}".format(i+1) + ": " + "{:10.2f}".format(self.Freq[i]))
        """ Compute and then print the TED. """
        # print(self.L.shape)
        self.TED = np.multiply(self.L,np.transpose(LA.pinv(self.L)))*100
        np.set_printoptions(suppress=True)
        print('Total Energy Distribution:')
        print(self.TED)


        tableOutput = "{:>26s}".format("Frequency #: ")
        for i in range(len(self.Freq)):
            tableOutput += "{:8d}".format(i+1)
        tableOutput += '\n'
        tableOutput += "--------------------------"
        for i in range(len(self.Freq)):
            tableOutput += "--------"
        tableOutput += '\n'
        tableOutput += "{:>26s}".format("Frequency: ")
        for i in range(len(self.Freq)):
            tableOutput += " " + "{:7.1f}".format(self.Freq[i])
        tableOutput += '\n'
        tableOutput += "--------------------------"
        for i in range(len(self.Freq)):
            tableOutput += "--------"
        tableOutput += '\n'
        """ Modify table rows so that internal coord labels have a fixed length. """
        for i in range(len(self.TED)):
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
            for j in range(len(self.TED[0])):
                tableOutput += "{:8.1f}".format(self.TED[i][j])
            tableOutput += '\n'
        print(tableOutput)

    def project(self):
        L_T = np.transpose(self.L)
        delArray = []
        for i in range(len(L_T)):
            if np.dot(L_T[i],np.transpose(L_T[i])) < self.projTol:
                delArray.append(i)
        if len(delArray):
            L_T = np.delete(L_T,delArray,0)
        self.L = np.transpose(L_T)
        if len(delArray):
            self.eig_v = np.delete(self.eig_v,delArray,0)

