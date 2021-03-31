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

    If this is the case, then A simply becomes the following:

    A = B^T(BB^T)^-1.

    Where in fact, BB^T = G (non-mass-weighted)

"""

class TransDisp(object):
    def __init__(self,s_vectors,zmat,disp,eigs,Conv,dispTol):
        self.dispTol         = dispTol
        self.Conv            = Conv
        self.s_vectors       = s_vectors
        self.zmat            = zmat
        self.refCarts        = zmat.CartesiansFinal.copy()
        self.refCarts        = np.array(self.refCarts).astype(float)
        self.u               = np.identity(3*len(zmat.atomList))
        self.disp            = disp
        # self.rdisp           = rdisp
        # self.adisp           = adisp
        self.eigs            = eigs
        self.DispCart        = {}
        self.DispCart["ref"] = self.refCarts.copy()

    def run(self,delArray):
        self.delArray = delArray
        self.B = self.s_vectors.B
        """ Invert the L-matrix and then normalize the rows. """
        self.eig_inv = LA.pinv(self.eigs.copy())
        for i in range(len(self.eig_inv)):
            self.eig_inv[i] = self.eig_inv[i]/LA.norm(self.eig_inv[i])
        """
            Construct 'A', the commented lines may be useful for getting intensities later. 
        """
        # self.A = (self.u).dot(self.B.transpose())
        # self.A = (self.B).dot(self.A)
        self.A = (self.B).dot(self.B.transpose())
        self.A = LA.inv(self.A)
        self.A = (self.B.transpose()).dot(self.A)
        """ Similar to the first two commented lines, this will be necessary for intensities. """
        # self.A = (self.u).dot(self.A)
        """ This step modifies A to convert from normal coords to carts. """
        self.A = np.dot(self.A,LA.pinv(self.eig_inv)).round(decimals=12)
        """
            Next, we will have to determine our desired Normal mode internal coordinate displacements
        """
        
        """ Generate the internal coordinate displacement vector """
        self.n_disp = np.array([])
        N_int = len(self.zmat.bondIndices) + len(self.zmat.angleIndices) + len(self.zmat.torsionIndices)
        """ Need to use one of the dimensions of eigenvectors to generate disp-vector """
        for i in range(len(self.eigs[0])):
            self.n_disp = np.append(self.n_disp,self.disp)
        
        self.n_coord = self.INTC(self.refCarts,self.eig_inv,self.delArray)

        self.dispSym = np.zeros(len(self.n_disp)-len(self.delArray))

        for i in range(len(self.n_disp)-len(self.delArray)):
            disp = np.zeros(len(self.n_disp)-len(self.delArray))
            disp[i] = self.disp
            self.DispCart[str(i+1)+'_plus'] = self.CoordConvert(disp,self.n_coord.copy(),self.refCarts.copy(),50,1.0e-14)
            self.DispCart[str(i+1)+'_minus'] = self.CoordConvert(-disp,self.n_coord.copy(),self.refCarts.copy(),50,1.0e-14)
            norm1 = LA.norm(self.DispCart[str(i+1)+'_plus'] - self.DispCart["ref"])
            norm2 = LA.norm(self.DispCart[str(i+1)+'_minus'] - self.DispCart["ref"])
            normDiff = np.abs(norm1-norm2)
            if self.dispTol > normDiff:
                self.dispSym[i] = 1
        
        self.dispSym = self.dispSym.astype(int)


    def INTC(self,carts,eig_inv,delArray):
        intCoord = np.array([])
        for i in range(len(self.zmat.bondIndices)):
            x1 = np.array(carts[int(self.zmat.bondIndices[i][0])-1]).astype(float)
            x2 = np.array(carts[int(self.zmat.bondIndices[i][1])-1]).astype(float)
            intCoord = np.append(intCoord,self.calcBond(x1,x2))
        for i in range(len(self.zmat.angleIndices)):
            x1 = np.array(carts[int(self.zmat.angleIndices[i][0])-1]).astype(float)
            x2 = np.array(carts[int(self.zmat.angleIndices[i][1])-1]).astype(float)
            x3 = np.array(carts[int(self.zmat.angleIndices[i][2])-1]).astype(float)
            a = self.calcAngle(x1,x2,x3)
            if self.Conv:
                Condition1 = float(self.zmat.variableDictionaryFinal[self.zmat.angleVariables[i]]) > 180.
                Condition2 = float(self.zmat.variableDictionaryFinal[self.zmat.angleVariables[i]]) < 0.
                if Condition1 or Condition2:
                    a = 2*np.pi - a
            intCoord = np.append(intCoord,a)
        for i in range(len(self.zmat.torsionIndices)):
            x1 = np.array(carts[int(self.zmat.torsionIndices[i][0])-1]).astype(float)
            x2 = np.array(carts[int(self.zmat.torsionIndices[i][1])-1]).astype(float)
            x3 = np.array(carts[int(self.zmat.torsionIndices[i][2])-1]).astype(float)
            x4 = np.array(carts[int(self.zmat.torsionIndices[i][3])-1]).astype(float)
            t = self.calcTors(x1,x2,x3,x4)
            if self.Conv:
                """ Modify these conditions to check computed result against  """
                Condition1 = float(self.zmat.variableDictionaryFinal[self.zmat.torsionVariables[i]]) > 135. and t*180./np.pi < -135.
                Condition2 = float(self.zmat.variableDictionaryFinal[self.zmat.torsionVariables[i]]) < -135. and t*180./np.pi > 135.
                if Condition1:
                    t += 2*np.pi
                if Condition2:
                    t -= 2*np.pi
            intCoord = np.append(intCoord,t)
        if len(delArray):
            intCoord = np.delete(intCoord,delArray)
        intCoord = np.dot(eig_inv,intCoord)
        return intCoord

    def calcBond(self,x1,x2):
        r = LA.norm(x1-x2)
        return r

    def calcAngle(self,x1,x2,x3):
        a = np.dot(x1-x2,x3-x2)/(self.calcBond(x1,x2)*self.calcBond(x3,x2))
        a = np.arccos(a)
        return a

    def calcTors(self,x1,x2,x3,x4):
        e1 = (x2-x1)/self.calcBond(x1,x2)
        e2 = (x3-x2)/self.calcBond(x2,x3)
        e3 = (x4-x3)/self.calcBond(x3,x4)
        s  = np.dot(-e1,np.cross(-e2,e3))
        c  = np.dot(np.cross(-e1,e2),np.cross(-e2,e3))
        t = np.arctan2(s,c)
        return t

    def calcOutOfPlane(self,x1,x2,x3,x4):
        """ 
            This function will compute an out of plane angle between one bond and a plane
            formed by 3 other atoms. See page 58 of Molecular vibrations by Wilson, Cecius, and Cross
            for more info.
        """
        e1 = (x1-x2)/self.calcBond(x1,x2)
        e3 = (x3-x2)/self.calcBond(x3,x2)
        e4 = (x4-x2)/self.calcBond(x4,x2)
        phi = self.calcAngle(x3,x2,x4)
        theta = np.arcsin(np.dot(np.cross(e3,e4)/np.sin(phi),e1))
        return theta

    def CoordConvert(self,n_disp,n_coord,refCarts,max_iter,tolerance):
        newN = n_coord + n_disp
        newCarts = np.array(refCarts).astype(float)
        for i in range(max_iter):
            cartDisp = np.dot(self.A,n_disp)
            cartDispShaped = np.reshape(cartDisp,(-1,3))
            newCarts += cartDispShaped
            coordCheck = self.INTC(newCarts,self.eig_inv,self.delArray)
            n_disp = newN - coordCheck
            if LA.norm(n_disp) < tolerance:
                break
        return newCarts

