import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from scipy.linalg import fractional_matrix_power
from scipy.linalg import block_diag

class GFMethod(object):
    """
    This class is a versatile method to compute GF frequencies.

    TODO: Insert standard uncertainties of amu_elmass and HARTREE_WAVENUM
    """

    #def __init__(self, G, F, tol, proj_tol, zmat, ted, cma=None):
    def __init__(self, G, F, tol, proj_tol, zmat, ted, options, symm_obj, algo, cma=None):
        self.G = G
        self.F = F
        self.tol = tol
        self.proj_tol = proj_tol
        self.zmat = zmat
        self.AMU_ELMASS = 5.48579909065 * (10 ** (-4))
        self.HARTREE_WAVENUM = 219474.6313708
        self.ted = ted
        self.symm_obj = symm_obj
        self.algo = algo
        self.cma = cma
        self.options = options
    def block_GF(self):
        self.eigvals, self.eigvecs = [], []
        offset_h = 0
        for h, irrep in enumerate(self.symm_obj.symtext.irreps):
            if self.symm_obj.symtext.irreps[h].d == 1:
                fo = np.zeros((self.symm_obj.proj_irreps[h],self.symm_obj.proj_irreps[h]))
                for index in self.algo.indices_by_irrep[h]:
                    fc = self.F_O[index[0], index[1]]
                    fo[index[0] - offset_h, index[1] - offset_h] = self.F_O[index[0], index[1]]
                    fo[index[1] - offset_h, index[0] - offset_h] = self.F_O[index[1], index[0]]
                offset_h += self.symm_obj.proj_irreps[h]
                eig_v_h, L_p_h = LA.eigh(fo)
                self.eigvals.append(eig_v_h)
                self.eigvecs.append(L_p_h)

            elif self.symm_obj.symtext.irreps[h].d > 1:
                for p in range(0, self.symm_obj.symtext.irreps[h].d):
                    fo = np.zeros((self.symm_obj.proj_irreps[h],self.symm_obj.proj_irreps[h]))
                    offset_pf = p * self.symm_obj.proj_irreps[h]
                    for index in self.algo.indices_by_irrep[h]:
                        fo[index[0] - offset_h, index[1] - offset_h] = self.F_O[index[0], index[1]]
                        fo[index[1] - offset_h, index[0] - offset_h] = self.F_O[index[1], index[0]]
                    eig_v_h, L_p_h = LA.eigh(fo)
                    self.eigvals.append(eig_v_h)
                    self.eigvecs.append(L_p_h)
                offset_h += self.symm_obj.proj_irreps[h] * self.symm_obj.symtext.irreps[h].d
        self.eigvals = [item for sublist in self.eigvals for item in sublist]
        self.eig_v = np.asarray(self.eigvals)
        self.L_p = block_diag(*self.eigvecs)

    def run(self):
        # Construct the orthogonalizer
        self.G_O = fractional_matrix_power(self.G, 0.5)

        # Symmetrize F, diagonalize, then backtransform the eigenvectors.
        self.F_O = np.dot(np.dot(self.G_O, self.F), self.G_O)
        if self.options.symmetry:
            self.block_GF()
        else:
            self.eig_v, self.L_p = LA.eigh(self.F_O)
        #self.eig_v, self.L_p = LA.eigh(self.F_O)
        self.L_p[np.abs(self.L_p) < self.tol] = 0
        self.L = np.dot(self.G_O, self.L_p)
        self.L = np.real(self.L)
        L = np.absolute(self.L)
        L_p = np.real(self.L_p)
        self.L_p = L_p
        # Construct the normal mode overlap matrix. Will be useful for off-diagonal diagnostics.
        L = np.absolute(np.real(self.L))
        L_inv = LA.inv(L)
        S = np.dot(L.T, L)
        self.S = S

        self.L[np.abs(self.L) < self.tol] = 0
        # Compute the frequencies by the square root of the eigenvalues.
        self.freq = np.sqrt(self.eig_v, dtype=complex)
        # Filter for imaginary modes.
        for i in range(len(self.freq)):
            if np.real(self.freq[i]) > 0.0:
                self.freq[i] = np.real(self.freq[i])
            else:
                self.freq[i] = -np.imag(self.freq[i])

        # This seems to fix the constant warning I get about casting out imaginary values from the complex numbers,
        # however I will need to test it on transition states to see if it captures the imaginary frequencies.
        self.freq = self.freq.astype(float)

        # Convert from Hartrees to wavenumbers.
        self.freq *= self.HARTREE_WAVENUM
        for i in range(len(self.freq)):
            print(
                "frequency #"
                + "{:3d}".format(i + 1)
                + ": "
                + "{:10.2f}".format(self.freq[i])
            )
        # Compute and then print the TED.
        print("////////////////////////////////////////////")
        print("//{:^40s}//".format("Total Energy Distribution (TED)"))
        print("////////////////////////////////////////////")
        self.ted.run(self.L, self.freq, rect_print=False)
        self.ted_breakdown = self.ted.ted_breakdown
