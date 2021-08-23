import numpy as np

"""
The purpose of this class is to return a list of indicies by which the force constants of the CMA method 
will be computed. These indices will be determined by user input or by a scoring function which takes into
consideration the overlap of the normal coordinates and the difference in force constants for a particular normal
mode.
"""

class Algorithm(object):
    def __init__(self,eigs,S,initial_fc,options):
    
        self.eigs = eigs
        self.options = options
        self.S = S
        self.initial_fc = initial_fc

    def run(self):
        S = self.S
        initial_fc = self.initial_fc
         
        def loop(off_diag,lim):
                  
            indices = []
            Sum = 2
            for i in range(a):
                for j in range(i,i+off_diag):
                    if j > a -1:
                        break
                    else: 
                        if i == j: 
                           indices.append([i,j]) 
                        elif i != j:
                            if i > lim:
                                break
                            else:
                                indices.append([i,j]) 
                        Sum +=2
            return indices

        def coupling_diagnostic(a,overlap, initial_fc, tolerance):
            indices = []
            Diag = np.zeros((a,a))
            print("overlap")
            print(S)
            for x in range(a):
                for y in range(a):
                    if x ==y:
                        Diag[x,x] = 0
                    elif x !=y:
                        Diag[x,y] = (S[x,y] / (initial_fc[x] - initial_fc[y]) ) 

            Diag = np.absolute(Diag)
            print('Diagnostic Matrix')
            print(Diag) 
            for x in range(a):
                for y in range(a):
                    if x ==y:
                        indices.append([x,y]) 
                    if x != y:
                        if Diag[x,y] > tolerance:
                           indices.append([x,y]) 
                        else:
                            break
            return indices 
        tolerance = 5e-32
        a = self.eigs
        if self.options.mode_coupling_check:
            self.indices = coupling_diagnostic(a,S,initial_fc, tolerance)
        else: 
             if self.options.off_diag:
                 off_diag = self.options.off_diag_bands + 1
                 if self.options.off_diag_limit != False:
                     lim = self.options.off_diag_limit - 1 
                 else:
                     lim = a + 1
             else:
                 lim = a + 1
                 off_diag = 1
             self.indices = loop(off_diag,lim)
