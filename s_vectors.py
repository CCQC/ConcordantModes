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
    def __init__(self,zmat,options,variableDict):
        self.s_2center_dict    = {}
        self.s_3center_dict    = {}
        self.s_4center_dict    = {}
        self.bondIndices       = np.array(zmat.bondIndices).astype(np.int) 
        self.angleIndices      = np.array(zmat.angleIndices).astype(np.int) 
        self.torsionIndices    = np.array(zmat.torsionIndices).astype(np.int) 
        self.oopIndices        = np.array(zmat.oopIndices).astype(np.int) 
        self.options           = options
        self.variableDict      = variableDict

    def run(self,carts,B_proj):
        """ Initialize the cartesian coordinates """
        self.carts = carts
        """
            So, first things first I'll have to code up the 
            proper equations for the first order B-Tensors
            for:
            Bonds,
            Angles,
            Torsions,
            Out-of-plane bending,
            LinX Bends (Eventually),
            LinY Bends (Eventually)
        """
        # Toggle this to convert from bohr to angstrom
        # self.carts = 0.5291772085936*self.carts
        """
            First, bonds.
        """
        if len(self.bondIndices) > 0:
            for i in range(len(self.bondIndices)):
                self.s_2center_dict["B" + str(i+1)] = self.carts.copy()
                self.s_2center_dict["B" + str(i+1)] = 0*self.s_2center_dict["B" + str(i+1)]
                r = self.compute_r(self.carts,self.bondIndices[i][0]-1,self.bondIndices[i][1]-1)
                self.s_2center_dict["B" + str(i+1)][self.bondIndices[i][0]-1] = self.compute_e(self.carts,self.bondIndices[i][0]-1,self.bondIndices[i][1]-1,r)
                self.s_2center_dict["B" + str(i+1)][self.bondIndices[i][1]-1] = -self.s_2center_dict["B" + str(i+1)][self.bondIndices[i][0]-1]

        """
            Next, angles.
        """
        if len(self.angleIndices) > 0:
            for i in range(len(self.angleIndices)):
                self.s_3center_dict["A" + str(i+1)] = self.carts.copy()
                self.s_3center_dict["A" + str(i+1)] = 0*self.s_3center_dict["A" + str(i+1)]
                r_1 = self.compute_r(self.carts,self.angleIndices[i][0]-1,self.angleIndices[i][1]-1)
                r_2 = self.compute_r(self.carts,self.angleIndices[i][1]-1,self.angleIndices[i][2]-1)
                e_1 = self.compute_e(self.carts,self.angleIndices[i][0]-1,self.angleIndices[i][1]-1,r_1)
                e_2 = self.compute_e(self.carts,self.angleIndices[i][2]-1,self.angleIndices[i][1]-1,r_2)
                phi = self.compute_phi(e_1,e_2)
                self.s_3center_dict["A" + str(i+1)][self.angleIndices[i][0]-1] = self.compute_BEND(e_1,e_2,phi,r_1)
                self.s_3center_dict["A" + str(i+1)][self.angleIndices[i][2]-1] = self.compute_BEND(e_2,e_1,phi,r_2)
                self.s_3center_dict["A" + str(i+1)][self.angleIndices[i][1]-1] = -self.s_3center_dict["A" + str(i+1)][self.angleIndices[i][0]-1] - self.s_3center_dict["A" + str(i+1)][self.angleIndices[i][2]-1]
        
        """
            Next, torsions.
        """
        if len(self.torsionIndices) > 0:
            for i in range(len(self.torsionIndices)):
                self.s_4center_dict["D" + str(i+1)] = self.carts.copy()
                self.s_4center_dict["D" + str(i+1)] = 0*self.s_4center_dict["D" + str(i+1)]
                r_1   = self.compute_r(self.carts,self.torsionIndices[i][0]-1,self.torsionIndices[i][1]-1)
                r_2   = self.compute_r(self.carts,self.torsionIndices[i][1]-1,self.torsionIndices[i][2]-1)
                r_3   = self.compute_r(self.carts,self.torsionIndices[i][2]-1,self.torsionIndices[i][3]-1)
                e_1   = self.compute_e(self.carts,self.torsionIndices[i][0]-1,self.torsionIndices[i][1]-1,r_1)
                e_2   = self.compute_e(self.carts,self.torsionIndices[i][1]-1,self.torsionIndices[i][2]-1,r_2)
                e_3   = self.compute_e(self.carts,self.torsionIndices[i][2]-1,self.torsionIndices[i][3]-1,r_3)
                phi_1 = self.compute_phi(e_1,-e_2)
                phi_2 = self.compute_phi(e_2,-e_3)
                self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][0]-1] = self.compute_TORS1(e_1,-e_2,phi_1,r_1)
                self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][3]-1] = self.compute_TORS1(-e_3,e_2,phi_2,r_3)
                self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][1]-1] = self.compute_TORS2(e_1,-e_2,-e_3,phi_1,phi_2,r_1,r_2)
                self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][2]-1] = -self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][0]-1] \
                                                                                   -self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][1]-1] \
                                                                                   -self.s_4center_dict["D" + str(i+1)][self.torsionIndices[i][3]-1] 

        """
            Now, out of plane bending.
        """
        if len(self.oopIndices) > 0:
            for i in range(len(self.oopIndices)):
                self.s_4center_dict["O" + str(i+1)] = self.carts.copy()
                self.s_4center_dict["O" + str(i+1)] = 0*self.s_4center_dict["O" + str(i+1)]
                r_1   = self.compute_r(self.carts,self.oopIndices[i][0]-1,self.oopIndices[i][1]-1)
                r_2   = self.compute_r(self.carts,self.oopIndices[i][2]-1,self.oopIndices[i][1]-1)
                r_3   = self.compute_r(self.carts,self.oopIndices[i][3]-1,self.oopIndices[i][1]-1)
                e_1   = self.compute_e(self.carts,self.oopIndices[i][0]-1,self.oopIndices[i][1]-1,r_1)
                e_2   = self.compute_e(self.carts,self.oopIndices[i][2]-1,self.oopIndices[i][1]-1,r_2)
                e_3   = self.compute_e(self.carts,self.oopIndices[i][3]-1,self.oopIndices[i][1]-1,r_3)
                phi   = self.compute_phi(e_2,e_3)
                theta = self.variableDict["O"+str(i+1)]*np.pi/180.
                self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][0]-1] = self.compute_OOP1(e_1,e_2,e_3,r_1,theta,phi)
                self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][2]-1] = self.compute_OOP2(e_1,e_2,e_3,r_2,theta,phi)
                self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][3]-1] = self.compute_OOP2(-e_1,e_2,e_3,r_3,theta,phi)
                self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][1]-1] = -self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][0]-1] \
                                                                               -self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][2]-1] \
                                                                               -self.s_4center_dict["O" + str(i+1)][self.oopIndices[i][3]-1] 

        """
            The last step will be to concatenate all of the s-vectors into a singular B-tensor, in order of stretches, then bends, then torsions.
            Note: I am going to modify this to hold all 2-center, 3-center, and 4-center internal coordinates.
        """
        self.B = np.array([self.s_2center_dict['B1'].flatten()])
        """
            Append stretches
        """
        for i in range(len(self.s_2center_dict)-1):
            self.B = np.append(self.B,np.array([self.s_2center_dict['B'+str(i+2)].flatten()]),axis=0)
        """
            Append bends
        """
        for i in range(len(self.s_3center_dict)):
            self.B = np.append(self.B,np.array([self.s_3center_dict['A'+str(i+1)].flatten()]),axis=0)
        """
            Append torsions
        """
        for i in range(len(self.torsionIndices)):
            self.B = np.append(self.B,np.array([self.s_4center_dict['D'+str(i+1)].flatten()]),axis=0)
        """
            Append oop bends
        """
        for i in range(len(self.oopIndices)):
            self.B = np.append(self.B,np.array([self.s_4center_dict['O'+str(i+1)].flatten()]),axis=0)

        tol = 1e-10
        """ Experimental code to acquire natural internal coordinates from diagonalized BB^T Matrix """
        if self.options.coords.upper() != "ZMAT":
            Proj,eigs,_ = LA.svd(self.B)
            Proj[np.abs(Proj) < tol] = 0
            print("Proj singular values:")
            print(eigs)
            if B_proj:
                projArray = np.array(np.where(np.abs(eigs) > tol))
                self.Proj = Proj.T[:len(projArray[0])]
                self.Proj = self.Proj.T
                # for i in range(len(self.Proj.T)):
                    # self.Proj.T[i][np.abs(self.Proj.T[i]) < np.max(np.abs(self.Proj.T[i]))*proj_tol] = 0
        else:
            self.Proj = np.eye(len(self.B))
        
        """
            self.Proj may be used to transfrom from full set of internal coords to symmetrized internal coords.
            self.Proj.T may be used to transform from the symmetrized set to the full set of internal coords.
        """

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

    """
        See Wilson, Decius, and Cross' "Molecular Vibrations" page 60 for the OOP1 and OOP2 formulae
    """
    def compute_OOP1(self,e_1,e_2,e_3,r,theta,phi):
        s = (np.cross(e_2,e_3)/(np.cos(theta)*np.sin(phi)) - np.tan(theta)*e_1)/r
        return s
    
    def compute_OOP2(self,e_1,e_2,e_3,r,theta,phi):
        s = (np.cross(e_3,e_1)/(np.cos(theta)*np.sin(phi)) \
                - (np.tan(theta)*(e_2-np.cos(phi)*e_3))/(np.sin(phi)**2))/r
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

