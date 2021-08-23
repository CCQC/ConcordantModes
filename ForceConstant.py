import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA

"""
    This script will calculate the force constants of the CMA normal
    modes using numerical differentiation.

    [f(x+h) - 2f(x) + f(x-h)] / h^2 = 
    [E_(i_plus) - 2*E_(ref) + E_(i_minus)] / disp_size^2

    [f(x+h,y+k) - f(x+h,y) - f(x,y+k) + 2f(x,y) - f(x-h,y) - f(x,y-k) + f(x-h,y-k) / 2*h^2
    [E_(ij_plus) - E_(i_plus,j) - E_(i,j_plus) + E_(ref) - E_(i_minus,j) - E_(i,j_minus) + E_(i_minus,j_minus) / 2*disp_size^2
"""

class ForceConstant(object):
    def __init__(self,disp,p_en_array,m_en_array,ref_en,options,indices):
        self.options = options
        self.disp = disp
        self.p_en_array = p_en_array
        self.m_en_array = m_en_array
        self.ref_en = ref_en
        self.indices = indices
    def run(self):
    
        """Functions for computing the diagonal and off-diagonal force constants """
        def diag_fc(e_p, e_m, e_r, disp):
            fc = (e_p - 2*e_r + e_m)/(disp**2) 
            return fc

        def off_diag_fc(e_pp, e_pi, e_pj, e_mi, e_mj, e_mm, e_r, disp):
            fc =  ((e_pp - e_pi - e_pj + 2*e_r - e_mi - e_mj + e_mm)/(2*(disp**2)))
            return fc

        indices = self.indices  
        p_en_array = self.p_en_array 
        m_en_array = self.m_en_array  
        disp = self.disp
        e_r = self.ref_en
        a = self.p_en_array.shape[0]
        self.FC = np.zeros((a,a))
        for index in indices:
            i,j = index[0], index[1] 
            e_pi, e_pj = p_en_array[i,i], p_en_array[j,j] 
            e_mi, e_mj = m_en_array[i,i], m_en_array[j,j] 
            e_pp, e_mm = p_en_array[i,j], m_en_array[i,j]  
            if i == j:
                self.FC[i,i] = diag_fc(e_pi, e_mi, e_r, disp.disp[i]) 
            elif i !=j:
                self.FC[i,j] = off_diag_fc(e_pp, e_pi, e_pj, e_mi, e_mj, e_mm, e_r, disp.disp[i])
        #Take advantage of FC[i,j] = FC[j,i]
        cf = np.triu_indices(a,1)
        il = (cf[1],cf[0])
        self.FC[il] = self.FC[cf]
     
 
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

#    def fcCalc(self,e_p,e_m,e_r,h):
#        fc = (e_p - 2*e_r + e_m)/(h**2)
#        return fc
#
#    def run(self):
#        off_diag = self.options.off_diag
#        p_en_array = self.p_en_array
#        m_en_array = self.m_en_array
