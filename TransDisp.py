import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA
from ConcordantModes.s_vectors import s_vectors


"""
    The purpose of this script is to perform an iterative transformation 
    between internal coordinate displacements and cartesian displacements, 
    using a first order, linear B-tensor transformation.
    
    It is necessary to construct an inverse B-tensor called 'A', where,
    
    A = uB^T(BuB^T)^-1.
    
    u can simply be the I_3N identity matrix, so long as frequency intensities 
    are not desired.

    If this is the case, then A simply becomes the following:

    A = B^T(BB^T)^-1.

    Where in fact, BB^T = G (non-mass-weighted)

"""

class TransDisp(object):
    def __init__(self,s_vectors,zmat,disp,eigs,Conv,dispTol,TED,options,
            GF=None):
        """
            As is generally the case for Programmers, GF = None by default.
        """
        self.dispTol         = dispTol
        self.Conv            = Conv
        self.s_vectors       = s_vectors
        self.zmat            = zmat
        self.refCarts        = zmat.CartesiansFinal.copy()
        self.refCarts        = np.array(self.refCarts).astype(float)
        self.u               = np.identity(3*len(zmat.atomList))
        self.disp            = disp
        self.TED             = TED
        self.eigs            = eigs
        self.options         = options
        self.DispCart        = {}
        self.DispCart["ref"] = self.refCarts.copy()
        if GF:
            self.GF = GF

    def run(self):
        self.B = self.s_vectors.B.copy() # (redundant internals (s) x cartesians (3N))
        """ Invert the L-matrix and then normalize the rows. """
        projTol = 1.0e-3
        self.eig_inv = LA.inv(self.eigs) # (Normal modes (Q) x Sym internals (S) )
        for i in range(len(self.eig_inv)):
            self.eig_inv[i] = self.eig_inv[i]/LA.norm(self.eig_inv[i])
            self.eig_inv[i][np.abs(self.eig_inv[i]) < \
                    np.max(np.abs(self.eig_inv[i]))*projTol] = 0
        
        """
            Compute the A-matrix to convert from normal coordinates to cartesians
        """
        self.A = self.ComputeA(self.s_vectors.B.copy(),self.TED.Proj,
                               self.eig_inv,self.zmat.massWeight)

        """ Generate the normal coordinate values for the reference structure """
        self.n_coord = self.INTC(self.refCarts,self.eig_inv,self.TED.Proj)
        
        print("Normal Coordinate Values:")
        for i in range(len(self.n_coord)):
            print("Normal Coordinate #{:<4n}: {: 3.3f}".format(
                i+1,self.n_coord[i]))
        
        """
            Next, we will have to specify our Normal mode internal coordinate 
            displacement sizes
        """
        self.dispSym = np.zeros(len(self.eigs))
        
        # L = inv(self.eig_inv)

        Disp = self.disp
        self.disp = []
        for i in range(len(self.n_coord)):
            self.disp.append(Disp)
        
        """
            This code will be useful for submitting the computation with 
            reduced displacements
        """
        # if self.options.reducedDisp:
            # self.Freq = inv(np.diag(self.GF.Freq.copy()))
            # redDisp = np.dot(self.disp,self.Freq)
            # redDisp = redDisp / min(redDisp)
            # redDisp = redDisp*self.disp
            # self.disp = redDisp

        # raise RuntimeError

        """ Now we actually generate the displacements """
        for i in range(len(self.eigs.T)):
            disp = np.zeros(len(self.eigs.T))
            disp[i] = self.disp[i]
            print("Disp #" + str(i+1))
            print(disp)
            self.DispCart[str(i+1)+'_plus'] = \
                self.CoordConvert(disp,self.n_coord.copy(),
                                  self.refCarts.copy(),50,1.0e-11,
                                  self.A.copy(),self.options.tightDisp,
                                  self.zmat,self.options)
            self.DispCart[str(i+1)+'_minus'] = \
                self.CoordConvert(-disp,self.n_coord.copy(),
                                  self.refCarts.copy(),50,1.0e-11,
                                  self.A.copy(),self.options.tightDisp,
                                  self.zmat,self.options)
        

    def INTC(self,carts,eig_inv,Proj):
        """
            This is a function that computes all currently implemented and 
            specified internal coordinates from the desired cartesian 
            coordinates. The internal coordinate values will be appended 
            together in the order of the for-loops below. All other methods 
            in this package which use the INTC generated variables MUST
            abide by this ordering.
        """
        # tol = 1.0e-3
        intCoord = np.array([])
        for i in range(len(self.zmat.bondIndices)):
            x1 = np.array(
                    carts[int(self.zmat.bondIndices[i][0])-1]).astype(float)
            x2 = np.array(
                    carts[int(self.zmat.bondIndices[i][1])-1]).astype(float)
            intCoord = np.append(intCoord,self.calcBond(x1,x2))
        for i in range(len(self.zmat.angleIndices)):
            x1 = np.array(
                    carts[int(self.zmat.angleIndices[i][0])-1]).astype(float)
            x2 = np.array(
                    carts[int(self.zmat.angleIndices[i][1])-1]).astype(float)
            x3 = np.array(
                    carts[int(self.zmat.angleIndices[i][2])-1]).astype(float)
            a = self.calcAngle(x1,x2,x3)
            if self.Conv:
                Condition1 = float(
                        self.zmat.variableDictionaryFinal[
                            self.zmat.angleVariables[i]]) > 180.
                Condition2 = float(
                        self.zmat.variableDictionaryFinal[
                            self.zmat.angleVariables[i]]) < 0.
                if Condition1 or Condition2:
                    a = 2*np.pi - a
            intCoord = np.append(intCoord,a)
        for i in range(len(self.zmat.torsionIndices)):
            x1 = np.array(
                    carts[int(self.zmat.torsionIndices[i][0])-1]).astype(float)
            x2 = np.array(
                    carts[int(self.zmat.torsionIndices[i][1])-1]).astype(float)
            x3 = np.array(
                    carts[int(self.zmat.torsionIndices[i][2])-1]).astype(float)
            x4 = np.array(
                    carts[int(self.zmat.torsionIndices[i][3])-1]).astype(float)
            t = self.calcTors(x1,x2,x3,x4)
            if self.Conv:
                Condition1 = float(
                        self.zmat.variableDictionaryFinal[
                            self.zmat.torsionVariables[i]]) > 135. and (
                                    t*180./np.pi < -135.)
                Condition2 = float(
                        self.zmat.variableDictionaryFinal[
                            self.zmat.torsionVariables[i]]) < -135. and (
                                    t*180./np.pi > 135.)
                if Condition1:
                    t += 2*np.pi
                if Condition2:
                    t -= 2*np.pi
            intCoord = np.append(intCoord,t)
        for i in range(len(self.zmat.oopIndices)):
            x1 = np.array(
                    carts[int(self.zmat.oopIndices[i][0])-1]).astype(float)
            x2 = np.array(
                    carts[int(self.zmat.oopIndices[i][1])-1]).astype(float)
            x3 = np.array(
                    carts[int(self.zmat.oopIndices[i][2])-1]).astype(float)
            x4 = np.array(
                    carts[int(self.zmat.oopIndices[i][3])-1]).astype(float)
            o = self.calcOOP(x1,x2,x3,x4)
            if self.Conv:
                Condition1 = float(
                        self.zmat.variableDictionaryFinal[
                            self.zmat.oopVariables[i]]) > 180.
                Condition2 = float(
                        self.zmat.variableDictionaryFinal[
                            self.zmat.oopVariables[i]]) < -180.
                if Condition1:
                    o = o - 2*np.pi
                if Condition2:
                    o = o + 2*np.pi
            intCoord = np.append(intCoord,o)
        for i in range(len(self.zmat.linIndices)):
            x1 = np.array(
                    carts[int(self.zmat.linIndices[i][0])-1]).astype(float)
            x2 = np.array(
                    carts[int(self.zmat.linIndices[i][1])-1]).astype(float)
            x3 = np.array(
                    carts[int(self.zmat.linIndices[i][2])-1]).astype(float)
            x4 = np.array(
                    carts[int(self.zmat.linIndices[i][3])-1]).astype(float)
            l = self.calcLIN(x1,x2,x3,x4)
            intCoord = np.append(intCoord,l)
        intCoord = np.dot(Proj.T,intCoord)
        # intCoord[np.abs(intCoord) < tol*np.max(np.abs(intCoord))] = 0
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

    def calcOOP(self,x1,x2,x3,x4):
        """ 
            This function will compute an out of plane angle between one bond 
            and a plane formed by 3 other atoms. See page 58 of Molecular 
            vibrations by Wilson, Decius, and Cross for more info. However, 
            the indices have been rearranged as follows,
            1-->1
            2-->3
            3-->4
            4-->2.
        """
        e1 = (x1-x2)/self.calcBond(x1,x2)
        e3 = (x3-x2)/self.calcBond(x3,x2)
        e4 = (x4-x2)/self.calcBond(x4,x2)
        phi = self.calcAngle(x3,x2,x4)
        theta = np.arcsin(np.dot(np.cross(e3,e4)/np.sin(phi),e1))
        return theta
    
    def calcLIN(self,x1,x2,x3,x4):
        e1 = (x1-x2)/self.calcBond(x1,x2)
        e3 = (x3-x2)/self.calcBond(x3,x2)
        e4 = (x4-x2)/self.calcBond(x4,x2)
        theta = np.arcsin(np.dot(e4,np.cross(e3,e1)))
        return theta

    def CoordConvert(
            self,n_disp,n_coord,refCarts,max_iter,tolerance,A,tightDisp,zmat,
            options):
        tol = 1.0e-3
        convIter = 10
        newN = n_coord + n_disp
        newCarts = np.array(refCarts).astype(float)
        for i in range(max_iter):
            cartDisp = np.dot(n_disp,A)
            cartDispShaped = np.reshape(cartDisp,(-1,3))
            newCarts += cartDispShaped
            coordCheck = self.INTC(newCarts,self.eig_inv,self.TED.Proj)
            n_disp = newN - coordCheck
            n_disp[np.abs(n_disp) < tol*np.max(np.abs(n_disp))] = 0
            
            if tightDisp:
                sVec = s_vectors(zmat,options,zmat.variableDictionaryFinal)
                sVec.run(newCarts,False)
                A = self.ComputeA(sVec.B,self.TED.Proj,self.eig_inv,self.zmat.massWeight)
            if LA.norm(n_disp) < tolerance:
                break
        if LA.norm(n_disp) > tolerance:
            print("This displacement did not converge.")
            print("Norm:")
            print(LA.norm(n_disp))
            print("Tolerance:")
            print(tolerance)
        return newCarts

    def ComputeA(self,B,Proj,eig_inv,u):
        """
            Construct 'A', the commented lines may be useful for getting 
            intensities later. The BB^T product must be linearly independent 
            to work, so it will be projected into the symmetry adapted basis.
        """
        # BT = B.T # (3N x s)
        # A = u.dot(BT)
        # A = np.dot(B,BT) # (s x s)
        # A = LA.pinv(A) # (s x s)
        # A = BT.dot(A) # (3N x s)
        # A = np.dot(A,Proj).T # (S x 3N)
        
        L = inv(eig_inv)
        A = LA.pinv(B)
        A = np.dot(A,Proj) # (3N x S)
        A = A.T # (S x 3N)
       
        """ 
            Similar to the first two commented lines, this will be necessary 
            for intensities. 
        """
        # A = A.dot(np.sqrt(inv(u)))
        """ This step modifies A to convert from normal coords to carts. """
        A = np.dot(L.T,A) # (Q x 3N)

        return A
