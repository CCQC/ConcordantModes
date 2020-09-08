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
        # For now a test case (methanol) is used. But later I will feed these in to initialize.
        # They MUST be np arrays, however!
        # self.carts    = np.array([[-1.37093756,-0.02453364,0.00000403],
                                  # [1.30184940,0.12055381,0.00000330],
                                  # [-2.07199296,1.90293990,-0.00091896],
                                  # [-2.10348346,-0.97490875,1.67685231],
                                  # [-2.10340229,-0.97642099,-1.67602173],
                                  # [1.94110263,-1.57276879,-0.00001194]])
        # self.bondIndices    = np.array([[2,1],[3,1],[4,1],[5,1],[6,2]])
        # self.angleIndices   = np.array([[3,1,2],[4,1,2],[5,1,2],[6,2,1]])
        # self.torsionIndices = np.array([[4,1,2,3],[5,1,2,3],[6,2,1,3]])
        # self.carts    = np.array([[-0.0000000000,-1.4236653500,0.9926470200],
                                  # [0.0000000000,0.0000000000,-0.1250915700],
                                  # [-0.0000000000,1.4236653500,0.9926470200]])
        # self.bondIndices    = np.array([[2,1],[3,2]])
        # self.angleIndices   = np.array([[3,2,1]])
        # self.torsionIndices = np.array([])
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
        # An ad hoc check
        # print(self.carts)
        # print(self.bondIndices)
        # print(self.angleIndices)
        # print(self.torsionIndices)
        self.carts = 0.5291772085936*self.carts
        """
            First, bonds.
        """
        if len(self.bondIndices) > 0:
            print("Bond s-vectors:")
            for i in range(len(self.bondIndices)):
                self.s_STRE_dict["B" + str(i+1)] = self.carts.copy()
                self.s_STRE_dict["B" + str(i+1)] = 0*self.s_STRE_dict["B" + str(i+1)]
                r = self.compute_r(self.carts,self.bondIndices[i][0]-1,self.bondIndices[i][1]-1)
                self.s_STRE_dict["B" + str(i+1)][self.bondIndices[i][0]-1] = self.compute_e(self.carts,self.bondIndices[i][0]-1,self.bondIndices[i][1]-1,r)
                self.s_STRE_dict["B" + str(i+1)][self.bondIndices[i][1]-1] = -self.s_STRE_dict["B" + str(i+1)][self.bondIndices[i][0]-1]
                """
                    Print out s-vector matrix for internal coord
                """
                print("B" + str(i+1))
                print(self.s_STRE_dict["B" + str(i+1)])

        """
            Next, angles.
        """
        if len(self.angleIndices) > 0:
            print("Angle s-vectors:")
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
                print("A" + str(i+1))
                print(self.s_BEND_dict["A" + str(i+1)])
        
        """
            Finally, for now, torsions.
        """
        if len(self.torsionIndices) > 0:
            print("Torsion s-vectors:")
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
                print("D" + str(i+1))
                print(self.s_TORS_dict["D" + str(i+1)])

        # for key in self.s_STRE_dict:
            # print("STRE,BEND dots:")
            # for k in self.s_BEND_dict:
                # print(key)
                # print(k)
                # print(np.dot(self.s_STRE_dict[key],self.s_BEND_dict[k]).round(decimals=10))
            # print("STRE,TORS dots:")
            # for k in self.s_TORS_dict:
                # print(key)
                # print(k)
                # print(np.dot(self.s_STRE_dict[key],self.s_TORS_dict[k]).round(decimals=10))
        # for key in self.s_BEND_dict:
            # print("BEND,TORS dots:")
            # for k in self.s_TORS_dict:
                # print(key)
                # print(k)
                # print(np.dot(self.s_BEND_dict[key],self.s_TORS_dict[k]).round(decimals=10))

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

