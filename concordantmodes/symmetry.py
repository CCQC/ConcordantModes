import numpy as np
import scipy
from numpy import linalg as LA
from concordantmodes.masses import get_mass
import copy

np.printoptions(precision=3, suppress=True)
class Symmetry(object):
    
    def __init__(self, zmat, options):
        self.zmat = zmat
        self.options = options
    
    def run(self):
        if self.options.autosalcs:
            try:
                import molsym
                from molsym.salcs.internal_coordinates import InternalCoordinates
                from molsym.salcs.projection_op import ProjectionOp
            except ImportError:
                raise ImportError('MolSym library not detected in your environment')
            biggun = self.zmat.bond_variables + self.zmat.angle_variables + self.zmat.torsion_variables + self.zmat.oop_variables + self.zmat.linx_variables + self.zmat.liny_variables 
            schema = self.make_schema()
            self.ic_list = []
            for var in self.zmat.variables:
                new_coord = self.remove_one(self.zmat.index_dictionary[var])
                self.ic_list.append(new_coord)
            
            self.nbfxns = len(self.ic_list) 
            ics = []
            for i in range(len(self.ic_list)):
                ics.append([self.ic_list[i], biggun[i]])
            mol = molsym.Molecule.from_schema(schema)
            self.symtext = molsym.Symtext.from_molecule(mol)
            if self.options.subgroup:
                self.symtext = self.symtext.subgroup_symtext(self.options.subgroup)
            ICs = InternalCoordinates(self.symtext, ics)
            print(f"The Molecular Point Group is self.symtext.pg")
            print(f"If this was not expected, adjust the tolerances in the MolSym program or symmetrize your molecule.")
            self.salcs = ProjectionOp(self.symtext, ICs)
            self.process_salcs() 
        else:
            print("Make your own salcs via manual projection matrix")

    def partner_functions(self):
        new_salcs = []
        fxn_list = []
        new_salc_list = []
        for s, salc in enumerate(self.salcs.salc_sets):
            degen = self.symtext.irreps[s].d
            if salc is None:
                new_salcs.append(np.zeros((0, self.nbfxns)))
                fxn_list.append([])
            elif degen > 1:
                new = np.zeros((len(self.salcs.partner_function_sets_by_irrep[s]) * degen, self.nbfxns))
                for pf in range(0, degen):
                    for ss, Set in enumerate(self.salcs.partner_function_sets_by_irrep[s]):
                        cpi = pf * len(self.salcs.partner_function_sets_by_irrep[s]) + ss
                        function = self.salcs.partner_function_sets_by_irrep[s][ss][pf]
                        new[cpi, :] = self.salcs.salc_list[function].coeffs
                        new_salc_list.append(self.salcs.salc_list[function]) 
                fxn_list.append([1 for i in range(0, (len(salc) // degen))])
                new_salcs.append(new)    
            else:
                for ss, Set in enumerate(self.salcs.partner_function_sets_by_irrep[s]):
                    function = self.salcs.partner_function_sets_by_irrep[s][ss][0]
                    new_salc_list.append(self.salcs.salc_list[function])
                new_salcs.append(salc)
                fxn_list.append([1 for i in range(0, len(salc))])
        self.salcs.salc_list = new_salc_list
        self.salcs.salc_sets = new_salcs
        self.fxn_list = fxn_list
     
    def process_salcs(self):
        self.partner_functions() 
        self.irreplength = []
        if self.options.exploit_degen:
            for s, SET in enumerate(self.salcs.salc_sets):
                degen = self.symtext.irreps[s].d
                if SET.shape[0] == 0:
                    self.irreplength.append(0)
                else:
                    self.irreplength.append(SET.shape[0] // degen)
        else: 
            for SET in self.salcs.salc_sets:
                self.irreplength.append(SET.shape[0])

    def make_schema(self):
        qc_obj = {
            "symbols": self.zmat.atom_list,
            "mass" : [get_mass(x) for x in self.zmat.atom_list],
            "geometry" : self.zmat.cartesians_init.flatten().tolist(),
        }
        return qc_obj
    
    def remove_one(self, coordinate_index):
        new_indices = []
        for index in coordinate_index:
            new = int(index) - 1
            new_indices.append(new)
        return new_indices
    
    def package_salcs(self, salcs):
        sdict = dict.fromkeys(salcs.irreps)
        for s, salcs in enumerate(salcs.partner_functions[0]):
            s_irrep = salcs[0][1].irrep
            if len(salcs) > 1:
                if sdict[s_irrep] is None:
                    sdict[s_irrep] = [salc[1].coeffs for salc in salcs]
                else:
                    for s, salc in enumerate(salcs):
                        sdict[s_irrep][s] = np.vstack((sdict[s_irrep][s], salc[1].coeffs))
            elif len(salcs) == 1:
                for salc in salcs:
                    if sdict[s_irrep] is None:
                        sdict[s_irrep] = salc[1].coeffs
                    else:
                        sdict[s_irrep] = np.vstack((sdict[s_irrep], salc[1].coeffs))
        return sdict
    
    def make_proj_v2(self, s_vec):
        print("Symmetrize B-Matrix")
        self.b = []
        for salc in self.salcs.salc_sets:
            self.b.append(np.dot(salc, s_vec.B))
        print("Project Out Redundancies")
        if not self.options.qr_project:
            self.svd_project()
        else:
            self.qr_project()
    
    def qr_project(self):
        """
        Uses QR decomposition to detect linear dependencies in an array. If row 2 and 3 are linearly dependent, 
        this algorithm will remove index 3 as a default. However, this concept could be further used to project out 
        internal coordinates based on a scheme of coordinate hierarchy, molecular topology, etc...
        """
        Proj = []
        tol = 1e-5
        self.proj_irreps = []
        self.salcblocks = []
        for b, block in enumerate(self.b):
            if self.options.exploit_degen:
                if block.shape[0] == 0:
                    self.proj_irreps.append(0)
                elif self.symtext.irreps[b].d > 1:
                    Q, R = np.linalg.qr(block.T[:, :self.irreplength[b]])
                    rejects = []
                    for i in range(block[:self.irreplength[b], :].shape[0]):
                        if np.isclose(R[i,i], 0.0, atol=1e-5):
                            rejects.append(i)
                    offset = 0
                    for ir in range(0, self.symtext.irreps[b].d):
                        self.salcblocks.append(np.delete(self.salcs.salc_sets[b][offset:offset + self.irreplength[b]], rejects, 0))
                        offset += self.irreplength[b]
                    self.proj_irreps.append(self.salcblocks[b].shape[0])
                    print(f"These are the rejects to project out of the SALCs")
                    print(rejects)
                    print("projected salcs")
                    print(self.salcblocks[b])
                else:
                    Q, R = np.linalg.qr(block.T)
                    rejects = []
                    for i in range(block.shape[0]):
                        if np.isclose(R[i,i], 0.0, atol=1e-5):
                            rejects.append(i)
                    self.salcblocks.append(np.delete(self.salcs.salc_sets[b], rejects, 0))
                    print(f"These are the rejects to project out of the SALCs")
                    print(rejects)
                    print("projected salcs")
                    print(self.salcblocks[b])
                    self.proj_irreps.append(self.salcblocks[b].shape[0])
        self.bloc = np.row_stack((self.salcblocks)).T
    
    def svd_project(self):
        Proj = []
        tol = 1e-5
        self.proj_irreps = []
        self.salcblocks = []

        for b, block in enumerate(self.b):
            if self.options.exploit_degen:
                if block.shape[0] == 0:
                    self.proj_irreps.append(0)
                elif self.symtext.irreps[b].d > 1:
                    block_pf = block[:self.irreplength[b]]
                    proj, eigs, _ = LA.svd(block_pf)
                    proj[np.abs(proj) < tol] = 0
                    proj_array = np.array(np.where(np.abs(eigs) > tol))
                    newproj = proj.T[: len(proj_array[0])]
                    newproj = newproj.T
                    self.proj_irreps.append(newproj.shape[1])
                    offset = 0
                    for i in range(0, self.symtext.irreps[b].d):
                        self.salcblocks.append(np.dot(newproj.T, self.salcs.salc_sets[b][offset:offset + self.irreplength[b]]))
                        offset += self.irreplength[b]
                else:
                    proj, eigs, _ = LA.svd(block)
                    proj[np.abs(proj) < tol] = 0
                    proj_array = np.array(np.where(np.abs(eigs) > tol))
                    newproj = proj.T[: len(proj_array[0])]
                    newproj = newproj.T
                    self.proj_irreps.append(newproj.shape[1])
                    self.salcblocks.append(np.dot(newproj.T, self.salcs.salc_sets[b]))
        self.bloc = np.row_stack((self.salcblocks)).T

