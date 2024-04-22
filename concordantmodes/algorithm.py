import numpy as np
import os
import scipy
from scipy import stats


class Algorithm(object):
    """
    The purpose of this class is to return a list of indices by which the force constants of the CMA method
    will be computed. These indices will be determined by user input or by a scoring function which takes into
    consideration the overlap of the normal coordinates and the difference in force constants for a particular normal
    mode.
    """

    def __init__(self, a, options, symm_obj, diag = False):
        self.a = a
        self.options = options
        self.symm_obj = symm_obj
        self.diag = diag
    def run(self):
        if self.options.symmetry:
            if self.diag:
                self.loop_symm_diag()
            else:
                self.loop_symmetry()

        else:
            self.loop()

    def loop(self):
        self.indices = []
        for i in range():
            for j in range(i):
                self.indices.append([i, j])

    def loop_symm_diag(self):
        self.indices = []
        self.indices_by_irrep = []
        offset = 0
        for h, irrep in enumerate(self.symm_obj.symtext.irreps):
            irrep_indices = []
            if self.symm_obj.proj_irreps[h] == 0:
                self.indices_by_irrep.append([])
            elif self.symm_obj.symtext.irreps[h].d == 1:
                for i in range(offset, self.symm_obj.proj_irreps[h] + offset):
                    irrep_indices.append([i,i])
                self.indices_by_irrep.append(irrep_indices)
            else:
                for i in range(offset, self.symm_obj.proj_irreps[h] + offset):
                    irrep_indices.append([i,i])
                self.indices_by_irrep.append(irrep_indices)
            offset += self.symm_obj.proj_irreps[h] * self.symm_obj.symtext.irreps[h].d
        self.indices = [item for sublist in self.indices_by_irrep for item in sublist]
        

    def loop_symmetry(self):
        self.indices = []
        self.indices_by_irrep = []
        offset = 0
        for h, irrep in enumerate(self.symm_obj.symtext.irreps):
            irrep_indices = []
            if self.symm_obj.proj_irreps[h] == 0:
                self.indices_by_irrep.append([])
            elif self.symm_obj.symtext.irreps[h].d == 1:
                for i in range(offset, self.symm_obj.proj_irreps[h] + offset):
                    for j in range(i, self.symm_obj.proj_irreps[h] + offset):
                        irrep_indices.append([i,j])
                self.indices_by_irrep.append(irrep_indices)
            else:
                for i in range(offset, self.symm_obj.proj_irreps[h] + offset):
                    for j in range(i, self.symm_obj.proj_irreps[h] + offset):
                        irrep_indices.append([i,j])
                self.indices_by_irrep.append(irrep_indices)
            offset += self.symm_obj.proj_irreps[h] * self.symm_obj.symtext.irreps[h].d
        self.indices = [item for sublist in self.indices_by_irrep for item in sublist]
