import os
import shutil
import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA


class SVectors(object):
    """
    s-vectors: s_a^k = (B_ax^k,B_ay^k,B_az^k)
    Where a refers to the atom #, and k refers to the internal coordinate.
    So, all B-tensors can be contained within the s-vector
    for each atom in the molecular system of interest.
    """

    def __init__(self, zmat, options, variable_dict):
        self.s_2center_dict = {}
        self.s_3center_dict = {}
        self.s_4center_dict = {}
        self.bond_indices = np.array(zmat.bond_indices).astype(np.int)
        self.angle_indices = np.array(zmat.angle_indices).astype(np.int)
        self.torsion_indices = np.array(zmat.torsion_indices).astype(np.int)
        self.oop_indices = np.array(zmat.oop_indices).astype(np.int)
        self.lin_indices = np.array(zmat.lin_indices).astype(np.int)
        self.linx_indices = np.array(zmat.linx_indices).astype(np.int)
        self.liny_indices = np.array(zmat.liny_indices).astype(np.int)
        self.options = options
        self.variable_dict = variable_dict

    def run(self, carts, B_proj, proj=None):
        self.proj = proj
        # Initialize the cartesian coordinates
        self.carts = carts
        # So, first things first I'll have to code up the
        # proper equations for the first order B-Tensors
        # for:
        # Bonds,
        # Angles,
        # Torsions,
        # Out-of-plane bending,
        # Linear Bends
        # Toggle this to convert from bohr to angstrom
        # self.carts = 0.5291772085936*self.carts

        # First, bonds.
        if len(self.bond_indices) > 0:
            for i in range(len(self.bond_indices)):
                self.s_2center_dict["B" + str(i + 1)] = self.carts.copy()
                self.s_2center_dict["B" + str(i + 1)] = (
                    0 * self.s_2center_dict["B" + str(i + 1)]
                )
                r = self.compute_r(
                    self.carts, self.bond_indices[i][0] - 1, self.bond_indices[i][1] - 1
                )
                self.s_2center_dict["B" + str(i + 1)][
                    self.bond_indices[i][0] - 1
                ] = self.compute_e(
                    self.carts,
                    self.bond_indices[i][0] - 1,
                    self.bond_indices[i][1] - 1,
                    r,
                )
                self.s_2center_dict["B" + str(i + 1)][
                    self.bond_indices[i][1] - 1
                ] = -self.s_2center_dict["B" + str(i + 1)][self.bond_indices[i][0] - 1]

        # Next, angles.
        if len(self.angle_indices) > 0:
            for i in range(len(self.angle_indices)):
                self.s_3center_dict["A" + str(i + 1)] = self.carts.copy()
                self.s_3center_dict["A" + str(i + 1)] = (
                    0 * self.s_3center_dict["A" + str(i + 1)]
                )
                r_1 = self.compute_r(
                    self.carts,
                    self.angle_indices[i][0] - 1,
                    self.angle_indices[i][1] - 1,
                )
                r_2 = self.compute_r(
                    self.carts,
                    self.angle_indices[i][1] - 1,
                    self.angle_indices[i][2] - 1,
                )
                e_1 = self.compute_e(
                    self.carts,
                    self.angle_indices[i][0] - 1,
                    self.angle_indices[i][1] - 1,
                    r_1,
                )
                e_2 = self.compute_e(
                    self.carts,
                    self.angle_indices[i][2] - 1,
                    self.angle_indices[i][1] - 1,
                    r_2,
                )
                phi = self.compute_phi(e_1, e_2)
                a_1 = self.compute_BEND(e_1, e_2, phi, r_1)
                a_2 = self.compute_BEND(e_2, e_1, phi, r_2)
                self.s_3center_dict["A" + str(i + 1)][
                    self.angle_indices[i][0] - 1
                ] = a_1
                self.s_3center_dict["A" + str(i + 1)][
                    self.angle_indices[i][2] - 1
                ] = a_2
                self.s_3center_dict["A" + str(i + 1)][self.angle_indices[i][1] - 1] = (
                    -a_1 - a_2
                )

        # Next, torsions.
        if len(self.torsion_indices) > 0:
            for i in range(len(self.torsion_indices)):
                self.s_4center_dict["D" + str(i + 1)] = self.carts.copy()
                self.s_4center_dict["D" + str(i + 1)] = (
                    0 * self.s_4center_dict["D" + str(i + 1)]
                )
                r_1 = self.compute_r(
                    self.carts,
                    self.torsion_indices[i][0] - 1,
                    self.torsion_indices[i][1] - 1,
                )
                r_2 = self.compute_r(
                    self.carts,
                    self.torsion_indices[i][1] - 1,
                    self.torsion_indices[i][2] - 1,
                )
                r_3 = self.compute_r(
                    self.carts,
                    self.torsion_indices[i][2] - 1,
                    self.torsion_indices[i][3] - 1,
                )
                e_1 = self.compute_e(
                    self.carts,
                    self.torsion_indices[i][0] - 1,
                    self.torsion_indices[i][1] - 1,
                    r_1,
                )
                e_2 = self.compute_e(
                    self.carts,
                    self.torsion_indices[i][1] - 1,
                    self.torsion_indices[i][2] - 1,
                    r_2,
                )
                e_3 = self.compute_e(
                    self.carts,
                    self.torsion_indices[i][2] - 1,
                    self.torsion_indices[i][3] - 1,
                    r_3,
                )
                phi_1 = self.compute_phi(e_1, -e_2)
                phi_2 = self.compute_phi(e_2, -e_3)
                t_1 = self.compute_TORS1(e_1, -e_2, phi_1, r_1)
                t_4 = self.compute_TORS1(-e_3, e_2, phi_2, r_3)
                t_2 = self.compute_TORS2(e_1, -e_2, -e_3, phi_1, phi_2, r_1, r_2)
                self.s_4center_dict["D" + str(i + 1)][
                    self.torsion_indices[i][0] - 1
                ] = t_1
                self.s_4center_dict["D" + str(i + 1)][
                    self.torsion_indices[i][3] - 1
                ] = t_4
                self.s_4center_dict["D" + str(i + 1)][
                    self.torsion_indices[i][1] - 1
                ] = t_2
                self.s_4center_dict["D" + str(i + 1)][
                    self.torsion_indices[i][2] - 1
                ] = (-t_1 - t_2 - t_4)

        # Now, out of plane bending.
        if len(self.oop_indices) > 0:
            for i in range(len(self.oop_indices)):
                self.s_4center_dict["O" + str(i + 1)] = self.carts.copy()
                self.s_4center_dict["O" + str(i + 1)] = (
                    0 * self.s_4center_dict["O" + str(i + 1)]
                )
                r_1 = self.compute_r(
                    self.carts, self.oop_indices[i][0] - 1, self.oop_indices[i][1] - 1
                )
                r_2 = self.compute_r(
                    self.carts, self.oop_indices[i][2] - 1, self.oop_indices[i][1] - 1
                )
                r_3 = self.compute_r(
                    self.carts, self.oop_indices[i][3] - 1, self.oop_indices[i][1] - 1
                )
                e_1 = self.compute_e(
                    self.carts,
                    self.oop_indices[i][0] - 1,
                    self.oop_indices[i][1] - 1,
                    r_1,
                )
                e_2 = self.compute_e(
                    self.carts,
                    self.oop_indices[i][2] - 1,
                    self.oop_indices[i][1] - 1,
                    r_2,
                )
                e_3 = self.compute_e(
                    self.carts,
                    self.oop_indices[i][3] - 1,
                    self.oop_indices[i][1] - 1,
                    r_3,
                )
                phi = self.compute_phi(e_2, e_3)
                theta = self.variable_dict["O" + str(i + 1)] * np.pi / 180.0
                o_1 = self.compute_OOP1(e_1, e_2, e_3, r_1, theta, phi)
                o_3 = self.compute_OOP2(e_1, e_2, e_3, r_2, theta, phi)
                o_4 = self.compute_OOP2(-e_1, e_3, e_2, r_3, theta, phi)
                self.s_4center_dict["O" + str(i + 1)][self.oop_indices[i][0] - 1] = o_1
                self.s_4center_dict["O" + str(i + 1)][self.oop_indices[i][2] - 1] = o_3
                self.s_4center_dict["O" + str(i + 1)][self.oop_indices[i][3] - 1] = o_4
                self.s_4center_dict["O" + str(i + 1)][self.oop_indices[i][1] - 1] = (
                    -o_1 - o_3 - o_4
                )

        # Linear bending.
        if len(self.lin_indices) > 0:
            for i in range(len(self.lin_indices)):
                self.s_4center_dict["L" + str(i + 1)] = self.carts.copy()
                self.s_4center_dict["L" + str(i + 1)] = (
                    0 * self.s_4center_dict["L" + str(i + 1)]
                )
                r_1 = self.compute_r(
                    self.carts, self.lin_indices[i][0] - 1, self.lin_indices[i][1] - 1
                )
                r_2 = self.compute_r(
                    self.carts, self.lin_indices[i][2] - 1, self.lin_indices[i][1] - 1
                )
                r_3 = self.compute_r(
                    self.carts, self.lin_indices[i][3] - 1, self.lin_indices[i][1] - 1
                )
                e_1 = self.compute_e(
                    self.carts,
                    self.lin_indices[i][0] - 1,
                    self.lin_indices[i][1] - 1,
                    r_1,
                )
                e_2 = self.compute_e(
                    self.carts,
                    self.lin_indices[i][2] - 1,
                    self.lin_indices[i][1] - 1,
                    r_2,
                )
                e_3 = self.compute_e(
                    self.carts,
                    self.lin_indices[i][3] - 1,
                    self.lin_indices[i][1] - 1,
                    r_3,
                )
                theta = self.variable_dict["L" + str(i + 1)] * np.pi / 180.0
                l_1 = self.compute_LIN(e_1, e_2, e_3, r_1, theta)
                l_3 = self.compute_LIN(e_2, e_3, e_1, r_2, theta)
                l_4 = self.compute_LIN(e_3, e_1, e_2, r_3, theta)
                self.s_4center_dict["L" + str(i + 1)][self.lin_indices[i][0] - 1] = l_1
                self.s_4center_dict["L" + str(i + 1)][self.lin_indices[i][2] - 1] = l_3
                self.s_4center_dict["L" + str(i + 1)][self.lin_indices[i][3] - 1] = l_4
                self.s_4center_dict["L" + str(i + 1)][self.lin_indices[i][1] - 1] = (
                    -l_1 - l_3 - l_4
                )

        # LinX bending.
        if len(self.linx_indices) > 0:
            for i in range(len(self.linx_indices)):
                self.s_4center_dict["Lx" + str(i + 1)] = self.carts.copy()
                self.s_4center_dict["Lx" + str(i + 1)] = (
                    0 * self.s_4center_dict["Lx" + str(i + 1)]
                )
                r_1 = self.compute_r(
                    self.carts, self.linx_indices[i][1] - 1, self.linx_indices[i][0] - 1
                )
                r_2 = self.compute_r(
                    self.carts, self.linx_indices[i][2] - 1, self.linx_indices[i][1] - 1
                )
                r_3 = self.compute_r(
                    self.carts, self.linx_indices[i][3] - 1, self.linx_indices[i][2] - 1
                )
                e_1 = self.compute_e(
                    self.carts,
                    self.linx_indices[i][0] - 1,
                    self.linx_indices[i][1] - 1,
                    r_1,
                )
                e_2 = self.compute_e(
                    self.carts,
                    self.linx_indices[i][1] - 1,
                    self.linx_indices[i][2] - 1,
                    r_2,
                )
                e_3 = self.compute_e(
                    self.carts,
                    self.linx_indices[i][2] - 1,
                    self.linx_indices[i][3] - 1,
                    r_3,
                )
                ax = self.calc_alpha_x(e_1, e_2, e_3)
                phi_1 = self.compute_phi(-e_1, e_2)
                phi_2 = self.compute_phi(-e_2, e_3)
                lx_1 = self.compute_LINX1(e_1, e_2, e_3, r_1, phi_1, phi_2, ax)
                lx_2 = self.compute_LINX2(e_1, e_2, e_3, r_1, r_2, phi_1, phi_2, ax)
                lx_4 = self.compute_LINX4(e_1, e_2, e_3, r_3, phi_1, ax)
                self.s_4center_dict["Lx" + str(i + 1)][
                    self.linx_indices[i][0] - 1
                ] = lx_1
                self.s_4center_dict["Lx" + str(i + 1)][
                    self.linx_indices[i][1] - 1
                ] = lx_2
                self.s_4center_dict["Lx" + str(i + 1)][
                    self.linx_indices[i][3] - 1
                ] = lx_4
                self.s_4center_dict["Lx" + str(i + 1)][self.linx_indices[i][2] - 1] = (
                    -lx_1 - lx_2 - lx_4
                )

        # LinY bending.
        if len(self.liny_indices) > 0:
            for i in range(len(self.liny_indices)):
                self.s_4center_dict["Ly" + str(i + 1)] = self.carts.copy()
                self.s_4center_dict["Ly" + str(i + 1)] = (
                    0 * self.s_4center_dict["Ly" + str(i + 1)]
                )
                r_1 = self.compute_r(
                    self.carts, self.liny_indices[i][1] - 1, self.liny_indices[i][0] - 1
                )
                r_2 = self.compute_r(
                    self.carts, self.liny_indices[i][2] - 1, self.liny_indices[i][1] - 1
                )
                r_3 = self.compute_r(
                    self.carts, self.liny_indices[i][3] - 1, self.liny_indices[i][2] - 1
                )
                e_1 = self.compute_e(
                    self.carts,
                    self.liny_indices[i][0] - 1,
                    self.liny_indices[i][1] - 1,
                    r_1,
                )
                e_2 = self.compute_e(
                    self.carts,
                    self.liny_indices[i][1] - 1,
                    self.liny_indices[i][2] - 1,
                    r_2,
                )
                e_3 = self.compute_e(
                    self.carts,
                    self.liny_indices[i][2] - 1,
                    self.liny_indices[i][3] - 1,
                    r_3,
                )
                ay = self.calc_alpha_y(e_1, e_2, e_3)
                phi_1 = self.compute_phi(-e_1, e_2)
                phi_2 = self.compute_phi(-e_2, e_3)
                ly_1 = self.compute_LINY1(e_1, e_2, e_3, r_1, phi_1, ay)
                ly_2 = self.compute_LINY2(e_1, e_2, e_3, r_1, r_2, phi_1, ay)
                ly_4 = self.compute_LINY4(e_1, e_2, e_3, r_2, r_3, phi_1, ay)
                self.s_4center_dict["Ly" + str(i + 1)][
                    self.liny_indices[i][0] - 1
                ] = ly_1
                self.s_4center_dict["Ly" + str(i + 1)][
                    self.liny_indices[i][1] - 1
                ] = ly_2
                self.s_4center_dict["Ly" + str(i + 1)][
                    self.liny_indices[i][3] - 1
                ] = ly_4
                self.s_4center_dict["Ly" + str(i + 1)][self.liny_indices[i][2] - 1] = (
                    -ly_1 - ly_2 - ly_4
                )

        # The last step will be to concatenate all of the s-vectors into a singular B-tensor, in order of stretches, then bends, then torsions.
        # Note: I am going to modify this to hold all 2-center, 3-center, and 4-center internal coordinates.
        self.B = np.array([self.s_2center_dict["B1"].flatten()])
        # Append stretches
        for i in range(len(self.s_2center_dict) - 1):
            self.B = np.append(
                self.B,
                np.array([self.s_2center_dict["B" + str(i + 2)].flatten()]),
                axis=0,
            )
        # Append bends
        for i in range(len(self.s_3center_dict)):
            self.B = np.append(
                self.B,
                np.array([self.s_3center_dict["A" + str(i + 1)].flatten()]),
                axis=0,
            )
        # Append torsions
        for i in range(len(self.torsion_indices)):
            self.B = np.append(
                self.B,
                np.array([self.s_4center_dict["D" + str(i + 1)].flatten()]),
                axis=0,
            )
        # Append oop bends
        for i in range(len(self.oop_indices)):
            self.B = np.append(
                self.B,
                np.array([self.s_4center_dict["O" + str(i + 1)].flatten()]),
                axis=0,
            )
        # Append lin bends
        for i in range(len(self.lin_indices)):
            self.B = np.append(
                self.B,
                np.array([self.s_4center_dict["L" + str(i + 1)].flatten()]),
                axis=0,
            )
        # Append linx bends
        for i in range(len(self.linx_indices)):
            self.B = np.append(
                self.B,
                np.array([self.s_4center_dict["Lx" + str(i + 1)].flatten()]),
                axis=0,
            )
        # Append liny bends
        for i in range(len(self.liny_indices)):
            self.B = np.append(
                self.B,
                np.array([self.s_4center_dict["Ly" + str(i + 1)].flatten()]),
                axis=0,
            )
        # raise RuntimeError

        tol = 1e-10
        # Now we acquire a linearly independant set of internal coordinates from the diagonalized
        # BB^T Matrix
        if not self.options.man_proj:
            if self.options.coords.upper() != "ZMAT":
                proj, eigs, _ = LA.svd(self.B)
                proj[np.abs(proj) < tol] = 0
                print("proj singular values:")
                print(eigs)
                if B_proj:
                    proj_array = np.array(np.where(np.abs(eigs) > tol))
                    self.proj = proj.T[: len(proj_array[0])]
                    self.proj = self.proj.T
                    # for i in range(len(self.proj.T)):
                    # self.proj.T[i][np.abs(self.proj.T[i]) < np.max(np.abs(self.proj.T[i]))*proj_tol] = 0
            else:
                self.proj = np.eye(len(self.B))
        # self.proj may be used to transfrom from full set of internal
        # coords to symmetrized internal coords. self.proj.T may be used
        # to transform from the symmetrized set to the full set of internal
        # coords.

        # Beware! The projected B matrix cannot be psuedo inverted to form
        # the A-matrix. You lose information.

    def compute_STRE(self, bond_indices, carts, r):
        s = (carts[bond_indices[0] - 1] - carts[bond_indices[1] - 1]) / r
        return s

    def compute_BEND(self, e_1, e_2, phi, r):
        s = (e_1 * np.cos(phi) - e_2) / (r * np.sin(phi))
        return s

    def compute_TORS1(self, e_1, e_2, phi, r):
        s = np.cross(e_1, e_2) / (r * np.sin(phi) ** 2)
        return s

    def compute_TORS2(self, e_1, e_2, e_3, phi_1, phi_2, r_1, r_2):
        s = ((r_2 - r_1 * np.cos(phi_1)) / (r_1 * r_2 * np.sin(phi_1) ** 2)) * np.cross(
            e_2, e_1
        ) + (np.cos(phi_2) / (r_2 * np.sin(phi_2) ** 2)) * np.cross(-e_2, e_3)
        return s

    # See Wilson, Decius, and Cross' "Molecular Vibrations" page 60 for the
    # OOP1 and OOP2 formulae.

    # The indices from the textbook have been rearranged in this
    # implementation as follows,
    # 1-->1
    # 2-->3
    # 3-->4
    # 4-->2.

    def compute_OOP1(self, e_1, e_2, e_3, r, theta, phi):
        s = (
            np.cross(e_2, e_3) / (np.cos(theta) * np.sin(phi)) - np.tan(theta) * e_1
        ) / r
        return s

    def compute_OOP2(self, e_1, e_2, e_3, r, theta, phi):
        s = (
            np.cross(e_3, e_1) / (np.cos(theta) * np.sin(phi))
            - (np.tan(theta) * (e_2 - np.cos(phi) * e_3)) / (np.sin(phi) ** 2)
        ) / r
        return s

    def compute_LIN(self, e_1, e_2, e_3, r, theta):
        s = (np.cross(e_3, e_2) - (np.dot(e_1, np.cross(e_3, e_2))) * e_1) / (
            np.cos(theta) * r
        )
        return s

    def compute_LINX1(self, e_1, e_2, e_3, r, phi_1, phi_2, ax):
        a = -ax * e_1 / (r * np.sin(phi_1) ** 2)
        b = -e_2 * (ax * (np.tan(phi_1) ** -1) + np.cos(phi_2)) / (r * np.sin(phi_1))
        c = -e_3 / (r * np.sin(phi_1))
        s = a + b + c
        return s

    def compute_LINX2(self, e_1, e_2, e_3, r_1, r_2, phi_1, phi_2, ax):
        a = (
            e_1
            * (ax * (r_2 - r_1 * np.cos(phi_1)) - r_1 * np.sin(phi_1) * np.cos(phi_2))
            / (r_1 * r_2 * np.sin(phi_1) ** 2)
        )
        b = (
            -e_2
            * (
                ax * (np.tan(phi_1) ** -1) * (r_1 * np.cos(phi_1) - r_2)
                + np.cos(phi_2) * (2 * r_1 * np.cos(phi_1) - r_2)
            )
            / (r_1 * r_2 * np.sin(phi_1))
        )
        c = -e_3 * (r_1 * np.cos(phi_1) - r_2) / (r_1 * r_2 * np.sin(phi_1))
        s = a + b + c
        return s

    def compute_LINX4(self, e_1, e_2, e_3, r, phi, ax):
        a = e_1 / (r * np.sin(phi))
        b = e_2 * (np.tan(phi) ** -1) / r
        c = ax * e_3 / r
        s = a + b + c
        return s

    def compute_LINY1(self, e_1, e_2, e_3, r, phi, ay):
        a = -ay * (np.tan(phi) ** -1) * (e_1 * np.cos(phi) + e_2) / (r * np.sin(phi))
        b = np.cross(e_2, -e_3) / (r * np.sin(phi))
        c = -e_1 * ay / r
        s = a + b + c
        return s

    def compute_LINY2(self, e_1, e_2, e_3, r_1, r_2, phi, ay):
        a = (
            -np.cos(phi)
            * ay
            * (e_1 * (r_1 - r_2 * np.cos(phi)) - e_2 * (r_2 - r_1 * np.cos(phi)))
            / (r_1 * r_2 * np.sin(phi) ** 2)
        )
        b = (r_2 * np.cross(e_2, -e_3) + r_1 * np.cross(-e_3, e_1)) / (
            r_1 * r_2 * np.sin(phi)
        )
        c = ay * (e_1 * r_2 + e_2 * r_1) / (r_1 * r_2)
        s = a + b + c
        return s

    def compute_LINY4(self, e_1, e_2, e_3, r_1, r_2, phi, ay):
        s = ((np.sin(phi) ** -1) * np.cross(e_1, e_2) + e_3 * ay) / r_2
        return s

    def compute_e(self, carts, ind1, ind2, r):
        e = (carts[ind1] - carts[ind2]) / r
        return e

    def compute_r(self, carts, ind1, ind2):
        r = LA.norm(carts[ind1] - carts[ind2])
        return r

    def compute_phi(self, e_1, e_2):
        phi = np.arccos(np.dot(e_1, e_2))
        return phi

    def calc_alpha_x(self, e_1, e_2, e_3):
        # e_1 = (x1 - x2) / self.calc_bond(x1, x2)
        # e_2 = (x3 - x2) / self.calc_bond(x2, x3)
        # e_3 = (x4 - x3) / self.calc_bond(x3, x4)
        theta = self.compute_phi(e_1, e_2)
        s = np.dot(np.cross(e_1, e_2), np.cross(-e_2, e_3))
        ax = s / np.sin(theta)
        return ax

    def calc_alpha_y(self, e_1, e_2, e_3):
        # e_1 = (x1 - x2) / self.calc_bond(x1, x2)
        # e_2 = (x3 - x2) / self.calc_bond(x2, x3)
        # e_3 = (x4 - x3) / self.calc_bond(x3, x4)
        theta = self.compute_phi(e_1, e_2)
        s = np.dot(e_1, np.cross(-e_2, e_3))
        ay = s / np.sin(theta)
        return ay
