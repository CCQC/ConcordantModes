import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

"""
    This class may be used to convert a cartesian F-matrix to an internal F-Matrix,
    and vice versa. 
"""

class F_conv(object):
    def __init__(self, F, s_vec, zmat):
        self.F = F
        self.Coord = Coord
        self.s_vec = s_vec
        self.zmat = zmat

    def run(self):
        """
            First construct the transpose of the A matrix. This could potentially be passed into the TransDisp.py
            class, rather than recomputing the A matrix. Or perhaps I could make the A-matrix its own module, 
            so that in the future intensities could be computed.
        """
        G = np.dot(self.s_vec.B,self.s_vec.B.transpose())
        self.A_T = np.dot(inv(G),self.s_vec.B)
        self.F = np.einsum('ia,jb,ab->ij',self.A_T,self.A_T,self.F)
        self.F = self.F/(0.5291772085936**2)
        self.F *= 4.35974394
        # print('F-Matrix:')
        # print(self.F)
