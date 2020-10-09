import os
import shutil
import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

"""
    s-vectors: s_a^k = (B_ax^k,B_ay^k,B_az^k)
    So, all B-tensors can be contained within the s-vector
    for each atom in the molecular system of interest.
"""


class s_vectors(object):
    def __init__(self,zmat):
        self.s_STRE_dict    = {}
        self.s_BEND_dict    = {}
        self.s_TORS_dict    = {}
        self.carts          = np.array(zmat.Cartesians).astype(np.float) 
        self.bondIndices    = np.array(zmat.bondIndices).astype(np.int) 
        self.angleIndices   = np.array(zmat.angleIndices).astype(np.int) 
        self.torsionIndices = np.array(zmat.torsionIndices).astype(np.int) 

    def run(self):
        """
            So, first things first I'll have to code up the 
            proper equations for the first order B-Tensors
            for:
            Bonds,
            Angles,
            Torsions,
            then out-of-plane motions (Eventually)
        """
        # Toggle this to convert from bohr to angstrom
        # self.carts = 0.5291772085936*self.carts
        """
            First, bonds.
        """
        if len(self.bondIndices) > 0:
            for i in range(len(self.bondIndices)):
                self.s_STRE_dict["B" + str(i+1)] = self.carts.copy()
                self.s_STRE_dict["B" + str(i+1)] = 0*self.s_STRE_dict["B" + str(i+1)]
                r = self.compute_r(self.carts,self.bondIndices[i][0]-1,self.bondIndices[i][1]-1)
                self.s_STRE_dict["B" + str(i+1)][self.bondIndices[i][0]-1] = self.compute_e(self.carts,self.bondIndices[i][0]-1,self.bondIndices[i][1]-1,r)
                self.s_STRE_dict["B" + str(i+1)][self.bondIndices[i][1]-1] = -self.s_STRE_dict["B" + str(i+1)][self.bondIndices[i][0]-1]

        """
            Next, angles.
        """
        if len(self.angleIndices) > 0:
            for i in range(len(self.angleIndices)):
                self.s_BEND_dict["A" + str(i+1)] = self.carts.copy()
                self.s_BEND_dict["A" + str(i+1)] = 0*self.s_BEND_dict["A" + str(i+1)]
                r_1 = self.compute_r(self.carts,self.angleIndices[i][0]-1,self.angleIndices[i][1]-1)
                r_2 = self.compute_r(self.carts,self.angleIndices[i][1]-1,self.angleIndices[i][2]-1)
                e_1 = self.compute_e(self.carts,self.angleIndices[i][0]-1,self.angleIndices[i][1]-1,r_1)
                e_2 = self.compute_e(self.carts,self.angleIndices[i][2]-1,self.angleIndices[i][1]-1,r_2)
                phi = self.compute_phi(e_1,e_2)
                self.s_BEND_dict["A" + str(i+1)][self.angleIndices[i][0]-1] = self.compute_BEND(e_1,e_2,phi,r_1)
                self.s_BEND_dict["A" + str(i+1)][self.angleIndices[i][2]-1] = self.compute_BEND(e_2,e_1,phi,r_2)
                self.s_BEND_dict["A" + str(i+1)][self.angleIndices[i][1]-1] = -self.s_BEND_dict["A" + str(i+1)][self.angleIndices[i][0]-1] - self.s_BEND_dict["A" + str(i+1)][self.angleIndices[i][2]-1]
        
        """
            Finally, for now, torsions.
        """
        if len(self.torsionIndices) > 0:
            for i in range(len(self.torsionIndices)):
                self.s_TORS_dict["D" + str(i+1)] = self.carts.copy()
                self.s_TORS_dict["D" + str(i+1)] = 0*self.s_TORS_dict["D" + str(i+1)]
                r_1   = self.compute_r(self.carts,self.torsionIndices[i][0]-1,self.torsionIndices[i][1]-1)
                r_2   = self.compute_r(self.carts,self.torsionIndices[i][1]-1,self.torsionIndices[i][2]-1)
                r_3   = self.compute_r(self.carts,self.torsionIndices[i][2]-1,self.torsionIndices[i][3]-1)
                e_1   = self.compute_e(self.carts,self.torsionIndices[i][0]-1,self.torsionIndices[i][1]-1,r_1)
                e_2   = self.compute_e(self.carts,self.torsionIndices[i][1]-1,self.torsionIndices[i][2]-1,r_2)
                e_3   = self.compute_e(self.carts,self.torsionIndices[i][2]-1,self.torsionIndices[i][3]-1,r_3)
                phi_1 = self.compute_phi(e_1,-e_2)
                phi_2 = self.compute_phi(e_2,-e_3)
                self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][0]-1] = self.compute_TORS1(e_1,-e_2,phi_1,r_1)
                self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][3]-1] = self.compute_TORS1(-e_3,e_2,phi_2,r_3)
                self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][1]-1] = self.compute_TORS2(e_1,-e_2,-e_3,phi_1,phi_2,r_1,r_2)
                self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][2]-1] = -self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][0]-1] \
                                                                                          -self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][1]-1] \
                                                                                          -self.s_TORS_dict["D" + str(i+1)][self.torsionIndices[i][3]-1] 

        """
            The last step will be to concatenate all of the s-vectors into a singular B-tensor, in order of stretches, then bends, then torsions.
        """
        self.B = np.array([self.s_STRE_dict['B1'].flatten()])
        """
            Append stretches
        """
        for i in range(len(self.s_STRE_dict)-1):
            self.B = np.append(self.B,np.array([self.s_STRE_dict['B'+str(i+2)].flatten()]),axis=0)
        """
            Append bends
        """
        for i in range(len(self.s_BEND_dict)):
            self.B = np.append(self.B,np.array([self.s_BEND_dict['A'+str(i+1)].flatten()]),axis=0)
        """
            Append torsions
        """
        for i in range(len(self.s_TORS_dict)):
            self.B = np.append(self.B,np.array([self.s_TORS_dict['D'+str(i+1)].flatten()]),axis=0)

    def compute_STRE(self,bondIndices,carts,r):
        s = (carts[bondIndices[0]-1] - carts[bondIndices[1]-1])/r
        return s

    def compute_BEND(self,e_1,e_2,phi,r):
        s = (e_1*np.cos(phi)-e_2)/(r*np.sin(phi))
        return s

    def compute_TORS1(self,e_1,e_2,phi,r):
        s = np.cross(e_1,e_2)/(r*np.sin(phi)**2)
        return s

    def compute_TORS2(self,e_1,e_2,e_3,phi_1,phi_2,r_1,r_2):
        s = ((r_2-r_1*np.cos(phi_1))/(r_1*r_2*np.sin(phi_1)**2))*np.cross(e_2,e_1) + \
            (np.cos(phi_2)/(r_2*np.sin(phi_2)**2))*np.cross(-e_2,e_3)
        return s

    def compute_e(self,carts,ind1,ind2,r):
        e = (carts[ind1] - carts[ind2])/r
        return e

    def compute_r(self,carts,ind1,ind2):
        r = LA.norm(carts[ind1] - carts[ind2])
        return r

    def compute_phi(self,e_1,e_2):
        phi = np.arccos(np.dot(e_1,e_2))
        return phi

