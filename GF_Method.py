import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power

"""
    This class will be used as a versatile method to compute GF frequencies!
"""

class GF_Method(object):
    def __init__(self, G, F, tol, projTol, zmat, TED):
        self.G               = G
        self.F               = F
        self.tol             = tol
        self.projTol         = projTol
        self.zmat            = zmat
        self.amu_elMass      = 5.48579909065*(10**(-4))
        self.Hartree_wavenum = 219474.6313708
        self.TED             = TED

    def run(self):
        """ Construct the orthogonalizer """
        self.G_O = fractional_matrix_power(self.G,0.5)
        """ Symmetrize F, diagonalize, then backtransform the eigenvectors. """
        self.F_O = np.dot(np.dot(self.G_O,self.F),self.G_O)
        self.eig_v,self.L_p = LA.eigh(self.F_O)
        self.L_p[np.abs(self.L_p) < self.tol] = 0
        self.L = np.dot(self.G_O,self.L_p)
        
        """ Construct the normal mode overlap matrix """
        L = np.absolute(np.real(self.L))
        L_inv = LA.inv(L) 
        L_t = L.T
        S = np.dot(L_t,L)
        self.S = S
       
        self.L[np.abs(self.L) < self.tol] = 0
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
        for i in range(len(self.Freq)):
            print("Frequency #" + "{:3d}".format(i+1) + ": " 
                  + "{:10.2f}".format(self.Freq[i]))
        """ Compute and then print the TED. """
        # np.set_printoptions(suppress=True)
        print('////////////////////////////////////////////')
        print("//{:^40s}//".format('Total Energy Distribution (TED)'))
        print('////////////////////////////////////////////')
        self.TED.run(self.L,self.Freq,rectPrint=False)
