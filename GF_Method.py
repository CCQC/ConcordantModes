import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power

"""
    This class will be used as a versatile method to compute GF frequencies!
"""

class GF_Method(object):
    def __init__(self, G, F, tol):
        self.G               = G
        self.F               = F
        self.tol             = tol
        self.amu_elMass      = 5.48579909065*(10**(-4))
        self.Hartree_wavenum = 219474.6313708

    def run(self):
        """ Construct the orthogonalizer """
        self.G_O = fractional_matrix_power(self.G,0.5)
        """ Symmetrize F, diagonalize, then backtransform the eigenvectors. """
        self.F_O = np.dot(np.dot(self.G_O,self.F),self.G_O)
        self.eig_v, self.L_p = LA.eigh(self.F_O)
        self.L_p[np.abs(self.L_p) < self.tol] = 0
        self.L = np.dot(self.G_O,self.L_p)
        self.L[np.abs(self.L**2) < self.tol] = 0
        self.L = np.real(self.L)
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
        print('Frequencies:')
        print(self.Freq)
        """ Compute and then print the TED. """
        self.TED = np.multiply(self.L,np.transpose(inv(self.L)))*100
        np.set_printoptions(suppress=True)
        print('Total Energy Distribution:')
        print(self.TED)
