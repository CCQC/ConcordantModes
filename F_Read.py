import numpy as np
import re
from numpy.linalg import inv
from numpy import linalg as LA

"""
    This class may be used to read in force constant information from the standard
    3-column formatted file. Idk what that format is called, but it's pretty standard.
"""

class F_Read(object):
    def __init__(self, FileName):
        self.FileName = FileName
        self.FC_regex = re.compile(r'(-?\d+\.\d+)')

    def run(self):
        with open(self.FileName,'r') as file:
            fc = file.read()

        self.FC_mat = re.findall(self.FC_regex,fc)
        self.FC_mat = np.reshape(self.FC_mat,(int(np.sqrt(len(self.FC_mat))),-1)).astype(float)
