import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

"""
    This class will be used to compute the G-matrix of the GF-Matrix Method.
    It should be constructed to give the flexibility of transformation via the L-matrix.
"""

class G_Matrix(object):
    def __init__(self, zmat, s_vectors):
        self.zmat      = zmat
        self.s_vectors = s_vectors

    def run(self,delArray):
        self.G = np.zeros((len(self.s_vectors.B)+len(delArray),len(self.s_vectors.B)+len(delArray)))
        
        b_len = len(self.zmat.bondIndices)
        a_len = len(self.zmat.angleIndices)
        t_len = len(self.zmat.torsionIndices)
        
        """
            These for-loops are constructed as such to take advantage of the inherent symmetry of the G-Matrix, 
            as well as how I have decided to store the s-vectors. It may work best to mass weight the s-vectors
            back in the s_vectors.py file so that they may be simply dotted in this step, however that is not
            necessary to the fundamental workings of the code for now.

        """
        for i in range(b_len):
            for j in range(b_len):
                overlap_indices = np.intersect1d(self.zmat.bondIndices[i],self.zmat.bondIndices[j]).astype(int)
                G_Mat_buff = self.compute_element(self.s_vectors.s_STRE_dict["B"+str(i+1)], \
                        self.s_vectors.s_STRE_dict["B"+str(j+1)],overlap_indices)
                self.G[i][j] = G_Mat_buff
                self.G[j][i] = G_Mat_buff
            for j in range(a_len):
                overlap_indices = np.intersect1d(self.zmat.bondIndices[i],self.zmat.angleIndices[j]).astype(int)
                G_Mat_buff = self.compute_element(self.s_vectors.s_STRE_dict["B"+str(i+1)], \
                        self.s_vectors.s_BEND_dict["A"+str(j+1)],overlap_indices)
                self.G[i][j+b_len] = G_Mat_buff
                self.G[j+b_len][i] = G_Mat_buff
            for j in range(t_len):
                overlap_indices = np.intersect1d(self.zmat.bondIndices[i],self.zmat.torsionIndices[j]).astype(int)
                G_Mat_buff = self.compute_element(self.s_vectors.s_STRE_dict["B"+str(i+1)], \
                        self.s_vectors.s_TORS_dict["D"+str(j+1)],overlap_indices)
                self.G[i][j+b_len+a_len] = G_Mat_buff
                self.G[j+b_len+a_len][i] = G_Mat_buff
        for i in range(a_len):
            for j in range(a_len):
                overlap_indices = np.intersect1d(self.zmat.angleIndices[i],self.zmat.angleIndices[j]).astype(int)
                G_Mat_buff = self.compute_element(self.s_vectors.s_BEND_dict["A"+str(i+1)], \
                        self.s_vectors.s_BEND_dict["A"+str(j+1)],overlap_indices)
                self.G[i+b_len][j+b_len] = G_Mat_buff
                self.G[j+b_len][i+b_len] = G_Mat_buff
            for j in range(t_len):
                overlap_indices = np.intersect1d(self.zmat.angleIndices[i],self.zmat.torsionIndices[j]).astype(int)
                G_Mat_buff = self.compute_element(self.s_vectors.s_BEND_dict["A"+str(i+1)], \
                        self.s_vectors.s_TORS_dict["D"+str(j+1)],overlap_indices)
                self.G[i+b_len][j+b_len+a_len] = G_Mat_buff
                self.G[j+b_len+a_len][i+b_len] = G_Mat_buff
        for i in range(t_len):
            for j in range(t_len):
                overlap_indices = np.intersect1d(self.zmat.torsionIndices[i],self.zmat.torsionIndices[j]).astype(int)
                G_Mat_buff = self.compute_element(self.s_vectors.s_TORS_dict["D"+str(i+1)], \
                        self.s_vectors.s_TORS_dict["D"+str(j+1)],overlap_indices)
                self.G[i+b_len+a_len][j+b_len+a_len] = G_Mat_buff
                self.G[j+b_len+a_len][i+b_len+a_len] = G_Mat_buff
        """
            Temper G
        """
        tol = 1e-12
        self.G[np.abs(self.G) < tol] = 0
    """
        Pretty straightforward function, computes the G-Matrix elements!
    """
    def compute_element(self,s_1,s_2,overlap_indices):
        G_Mat_Element = 0.
        for i in overlap_indices:
            if self.zmat.masses[i-1] > 0.1:
                G_Mat_Element += np.dot(s_1[i-1],s_2[i-1])/self.zmat.masses[i-1]
        return G_Mat_Element
