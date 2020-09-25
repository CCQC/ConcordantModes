import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power

"""
    This class will be used as a versatile method to compute GF frequencies!
"""

class GF_Method(object):
    # def __init__(self, L, G, F, tol, stage):
    def __init__(self, G, F, tol):
        # self.L   = L
        self.G   = G
        self.F   = F
        self.tol = tol
        # self.stage = stage

    def run(self):
        """ Construct the orthogonalizer """
        self.G_O = fractional_matrix_power(self.G,0.5)
        print(self.G_O)
        self.F_O = np.dot(np.dot(self.G_O,self.F),self.G_O)
        # self.F_O = np.dot(self.G,self.F)
        self.eig_v, self.L_p = LA.eigh(self.F_O)
        # self.eig_v, self.L = LA.eigh(self.F_O)
        self.L_p[np.abs(self.L_p) < self.tol] = 0
        self.L = np.dot(self.G_O,self.L_p)
        self.L[np.abs(self.L) < self.tol] = 0
        print('E-Vals (Normal Mode Force Constants):')
        print(self.eig_v)
        print('E-vecs:')
        print(self.L)
        self.Freq = np.sqrt(self.eig_v, dtype=np.complex)*1302.791
        for i in range(len(self.Freq)):
            if np.real(self.Freq[i]) > 0.:
                self.Freq[i] = np.real(self.Freq[i])
            else:
                self.Freq[i] = -np.imag(self.Freq[i])
        self.Freq = self.Freq.astype(float)
        print('Frequencies:')
        print(self.Freq)
        self.TED = np.multiply(self.L,np.transpose(inv(self.L)))*100
        np.set_printoptions(suppress=True)
        print('Total Energy Distribution:')
        print(self.TED)
