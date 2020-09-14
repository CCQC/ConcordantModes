import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA


"""
    The purpose of this script is to perform an iterative transformation between
    internal coordinate displacements and cartesian displacements, using a first
    order, linear B-tensor transformation.
    
    It is necessary to construct an inverse B-tensor called 'A', where,
    
    A = uB^T(BuB^T)^-1.
    
    u can simply be the I_3N identity matrix, so long as frequency intensities are
    not desired.
"""

class TransDisp(object):
    def __init__(self,s_vectors,zmat,rdisp,adisp,eigs):
        self.s_vectors = s_vectors
        self.B = self.s_vectors.B
        self.zmat = zmat
        self.refCarts = zmat.Cartesians.copy()
        self.u = np.identity(3*len(zmat.atomList))
        self.rdisp = rdisp
        self.adisp = adisp
        self.eigs = eigs
        self.DispCart = {}
        self.DispCart["ref"] = self.refCarts.copy()

    def run(self):
        """ Invert the L-matrix and then normalize the rows. """
        self.eig_inv = inv(self.eigs.copy())
        for i in range(len(self.eig_inv)):
            self.eig_inv[i] = self.eig_inv[i]/LA.norm(self.eig_inv[i])
        
        """
            Construct 'A' 
        """
        self.A = (self.u).dot(self.B.transpose())
        self.A = (self.B).dot(self.A)
        self.A = LA.inv(self.A)
        self.A = (self.B.transpose()).dot(self.A)
        self.A = (self.u).dot(self.A)
        """ This step modifies A to convert from normal coords to carts. """
        self.A = np.dot((self.A),inv(self.eig_inv)).round(decimals=12)
        """
            Next, we will have to determine our desired Normal mode internal coordinate displacements
        """
        
        """ Generate the internal coordinate displacement vector """
        self.s_disp = np.array([])
        self.s_coord = np.array([])
        for i in range(len(self.zmat.bondIndices)):
            self.s_disp = np.append(self.s_disp,self.rdisp)
            self.s_coord = np.append(self.s_coord,self.zmat.variableDictionary[self.zmat.bondVariables[i]])
        for i in range(len(self.zmat.angleIndices)):
            self.s_disp = np.append(self.s_disp,self.adisp)
            self.s_coord = np.append(self.s_coord,self.zmat.variableDictionary[self.zmat.angleVariables[i]])
        for i in range(len(self.zmat.torsionIndices)):
            self.s_disp = np.append(self.s_disp,self.adisp)
            self.s_coord = np.append(self.s_coord,self.zmat.variableDictionary[self.zmat.torsionVariables[i]])
        self.s_coord = self.s_coord.astype(float)
       
        """ Weight the internal coord displacement sizes """
        self.s_tensor = np.multiply(self.eig_inv,self.s_disp)

        """ Transform the internal coordinate displacements to normal mode internal coord disps. """
        self.n_disp = np.dot(self.s_tensor,self.s_disp)
        # self.n_disp = np.dot(self.eig_inv,self.s_disp)
        # self.n_coord = np.dot(self.eig_inv,self.s_coord)

        self.n_coord = self.INTC(self.refCarts)

        # for i in range(1):
        for i in range(len(self.n_disp)):
            disp = np.zeros(len(self.n_disp))
            disp[i] = self.rdisp
            # disp[i] = self.n_disp[i].copy()
            self.DispCart['+'+str(i)] = self.CoordConvert(disp,self.n_coord.copy(),self.refCarts.copy(),50,1.0e-10)
        
        raise RuntimeError

    def INTC(self,carts):
        intCoord = np.array([])
        for i in range(len(self.zmat.bondIndices)):
            x1 = np.array(carts[int(self.zmat.bondIndices[i][0])-1]).astype(float)
            x2 = np.array(carts[int(self.zmat.bondIndices[i][1])-1]).astype(float)
            intCoord = np.append(intCoord,self.calcBond(x1,x2))
        for i in range(len(self.zmat.angleIndices)):
            x1 = np.array(carts[int(self.zmat.angleIndices[i][0])-1]).astype(float)
            x2 = np.array(carts[int(self.zmat.angleIndices[i][1])-1]).astype(float)
            x3 = np.array(carts[int(self.zmat.angleIndices[i][2])-1]).astype(float)
            intCoord = np.append(intCoord,self.calcAngle(x1,x2,x3))
        for i in range(len(self.zmat.torsionIndices)):
            x1 = np.array(carts[int(self.zmat.torsionIndices[i][0])-1]).astype(float)
            x2 = np.array(carts[int(self.zmat.torsionIndices[i][1])-1]).astype(float)
            x3 = np.array(carts[int(self.zmat.torsionIndices[i][2])-1]).astype(float)
            x4 = np.array(carts[int(self.zmat.torsionIndices[i][3])-1]).astype(float)
            t = self.calcTors(x1,x2,x3,x4)
            Condition1 = float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) > np.pi/2 and float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) <= (3*np.pi)/2
            Condition2 = float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) < -np.pi/2 and float(self.zmat.variableDictionary[self.zmat.torsionVariables[i]]) >= -(3*np.pi)/2
            if Condition1 or Condition2:
                t += np.pi
            intCoord = np.append(intCoord,t)
        intCoord = np.dot(self.eig_inv,intCoord)
        return intCoord

    def calcBond(self,x1,x2):
        r = LA.norm(x1-x2)
        return r

    def calcAngle(self,x1,x2,x3):
        a = np.dot(x1-x2,x3-x2)/(self.calcBond(x1,x2)*self.calcBond(x3,x2))
        a = np.arccos(a)
        return a

    def calcTors(self,x1,x2,x3,x4):
        e1 = (x1-x2)/self.calcBond(x1,x2)
        e2 = (x3-x2)/self.calcBond(x2,x3)
        e3 = (x4-x3)/self.calcBond(x3,x4)
        t = np.dot(e1,np.cross(-e2,e3))/(np.dot(np.cross(e1,e2),np.cross(-e2,e3)))
        t = np.arctan(t)
        return t

    def CoordConvert(self,n_disp,n_coord,refCarts,max_iter,tolerance):
        newN = n_coord + n_disp
        # print(n_coord)
        # print(newN)
        newCarts = np.array(refCarts).astype(float)
        print(self.INTC(newCarts))
        print(n_coord)
        for i in range(max_iter):
            print(i)
            cartDisp = np.dot(self.A,n_disp)
            # print(n_disp)
            cartDispShaped = np.reshape(cartDisp,(-1,3))
            print(cartDispShaped)
            print(newCarts)
            newCarts += cartDispShaped
            print(newCarts)
            coordCheck = self.INTC(newCarts)
            # print(newN)
            # print(coordCheck)
            n_disp = newN - coordCheck
            if LA.norm(n_disp) < tolerance:
                break
        return newCarts

