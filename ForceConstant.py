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
    def __init__(self,disp,energiesDict):
        self.disp = disp
        self.energiesDict = energiesDict

    def run(self):
        self.FC = np.array([])
        e_r = self.energiesDict['ref']
        for i in range((len(self.energiesDict)-1)//2):
            e_p = self.energiesDict[str(i+1) + '_plus']
            e_m = self.energiesDict[str(i+1) + '_minus']
            self.FC = np.append(self.FC,self.fcCalc(e_p,e_m,e_r,self.disp.disp))

    def fcCalc(self,e_p,e_m,e_r,h):
        fc = (e_p - 2*e_r + e_m)/(h**2)
        return fc

