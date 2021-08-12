import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

"""
    This script will calculate the force constants of the CMA normal
    modes using numerical differentiation.

    [f(x+h) - 2f(x) + f(x-h)] / h^2 = 
    [E_(i_plus) - 2*E_(ref) + E_(i_minus)] / disp_size^2
"""

class ForceConstant(object):
    #def __init__(self,disp,energiesDict):
    def __init__(self,disp,p_en_array,m_en_array,ref_en,options):
        self.options = options
        self.disp = disp
        self.p_en_array = p_en_array
        self.m_en_array = m_en_array
        off_diag = self.options.off_diag        

        self.ref_en = ref_en
    
    def fcCalc(self,e_p,e_m,e_r,h):
        fc = (e_p - 2*e_r + e_m)/(h**2)
        return fc

    def run(self):
        off_diag = self.options.off_diag
        p_en_array = self.p_en_array
        m_en_array = self.m_en_array

        ref_en = self.ref_en
        e_r = ref_en
        dim = p_en_array.shape[0]
        self.FC = np.zeros((dim,dim))
        for i in range(dim):
            for j in range(dim):
                e_p = p_en_array[i,i] 
                e_m = m_en_array[i,i] 
                if i == j:
                    self.FC[i,i] = (e_p - 2*e_r + e_m)/(0.01**2) 
        for i in range(dim):
            for j in range(i, i + off_diag):
                if j > dim -1:
                    break
                else:
                    if i !=j:
                        e_pi = p_en_array[i,i] 
                        e_mi = m_en_array[i,i] 
                        e_pj = p_en_array[j,j] 
                        e_mj = m_en_array[j,j] 
                        e_pp = p_en_array[i,j]  
                        e_mm = m_en_array[i,j]  
                        self.FC[i,j] = ((e_pp - e_pi - e_pj + 2*e_r - e_mi - e_mj + e_mm)/(2*0.01*0.01))
        cf = np.triu_indices(dim,1)
        il = (cf[1],cf[0])
        self.FC[il] = self.FC[cf]
        print(self.FC)
 
#old diagonal only force constant code

#class ForceConstant(object):
#    def __init__(self,disp,energiesDict):
#        self.disp = disp
#        self.energiesDict = energiesDict
#
#    def run(self):
#        self.FC = np.array([])
#        e_r = self.energiesDict['ref']
#        for i in range((len(self.energiesDict)-1)//2):
#            e_p = self.energiesDict[str(i+1) + '_plus']
#            e_m = self.energiesDict[str(i+1) + '_minus']
#            self.FC = np.append(self.FC,self.fcCalc(e_p,e_m,e_r,self.disp.disp[i]))
#        # raise RuntimeError
#
#    def fcCalc(self,e_p,e_m,e_r,h):
#        fc = (e_p - 2*e_r + e_m)/(h**2)
#        return fc
