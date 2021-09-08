import time
import sys
import shutil
import os
import re
import subprocess
import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv
from scipy.linalg import fractional_matrix_power
from ConcordantModes.DirectoryTree import DirectoryTree
from ConcordantModes.F_Convert import F_conv
from ConcordantModes.F_Read import F_Read
from ConcordantModes.ForceConstant import ForceConstant
from ConcordantModes.G_Matrix import G_Matrix
from ConcordantModes.GF_Method import GF_Method
from ConcordantModes.Reap import Reap
from ConcordantModes.s_vectors import s_vectors
from ConcordantModes.TED import TED
from ConcordantModes.TransDisp import TransDisp
from ConcordantModes.vulcan_template import vulcan_template
from ConcordantModes.ZMAT import ZMAT
from ConcordantModes.int2cart import int2cart
from ConcordantModes.Algorithm import Algorithm


class ConcordantModes(object):
    def __init__(self, options):
        self.options = options
        """ 
            This constant is from:
            https://physics.nist.gov/cgi-bin/cuu/Value?hr 
            If this link dies, find the new link on NIST
            for the Hartree to Joule conversion and pop it in there.
            There is a standard uncertainty of 0.0000000000085 to this value.
        """
        self.mdyne_Hart = 4.3597447222071
        """
            Standard uncertainty of 0.00000000080
        """
        self.Bohr_Ang = 0.529177210903

    def run(self):
        t1 = time.time()

        rootdir = os.getcwd()

        """
            Parse the output to get all pertinent ZMAT info
        """
        self.zmat = ZMAT(self.options)
        self.zmat.run()
        if self.options.geomCheck:
            raise RuntimeError

        """
            Compute the initial s-vectors
        """
        s_vec = s_vectors(self.zmat, self.options, self.zmat.variableDictionaryInit)
        s_vec.run(self.zmat.CartesiansInit, True)

        self.TED = TED(s_vec.Proj, self.zmat)

        """
            Compute G-Matrix
        """
        g_mat = G_Matrix(self.zmat, s_vec, self.options)
        g_mat.run()

        G = g_mat.G.copy()

        """
            Read in FC matrix in cartesians, then convert to internals.
        
            I am adding a process to compute an initial hessian in internal coordinates.
        """
        initBool = False
        if os.path.exists(rootdir + "/fc.dat"):
            f_read = F_Read("fc.dat")
        elif os.path.exists(rootdir + "/FCMFINAL"):
            f_read = F_Read("FCMFINAL")
        else:
            initBool = True

            # First generate displacements in internal coordinates
            indices = np.triu_indices(len(s_vec.Proj.T))
            indices = np.array(indices).T
            eigsInit = np.eye(len(s_vec.Proj.T))

            initDisp = TransDisp(
                s_vec,
                self.zmat,
                self.options.disp,
                eigsInit,
                True,
                self.options.dispTol,
                self.TED,
                self.options,
                indices,
            )
            initDisp.run()
            progInit = self.options.programInit
            prognameInit = progInit.split("@")[0]

            DirObjInit = DirectoryTree(
                prognameInit,
                self.zmat,
                initDisp,
                self.options.cartInsertInit,
                initDisp.p_disp,
                initDisp.m_disp,
                self.options,
                indices,
                "templateInit.dat",
            )
            DirObjInit.run()
            os.chdir(rootdir + "/Disps")
            dispList = []
            for i in os.listdir(rootdir + "/Disps"):
                dispList.append(i)

            v_template = vulcan_template(
                self.options, len(dispList), prognameInit, progInit
            )
            out = v_template.run()
            with open("displacements.sh", "w") as file:
                file.write(out)

            pipe = subprocess.PIPE

            process = subprocess.run(
                "qsub displacements.sh", stdout=pipe, stderr=pipe, shell=True
            )
            self.outRegex = re.compile(r"Your\s*job\-array\s*(\d*)")
            self.job_id = int(re.search(self.outRegex, str(process.stdout)).group(1))

            self.jobFinRegex = re.compile(r"taskid")
            while True:
                qacct_proc = subprocess.run(
                    ["qacct", "-j", str(self.job_id)], stdout=pipe, stderr=pipe
                )
                qacct_string = str(qacct_proc.stdout)
                job_match = re.findall(self.jobFinRegex, qacct_string)
                if len(job_match) == len(dispList):
                    break
                time.sleep(10)

            output = str(process.stdout)
            error = str(process.stderr)

            ReapObjInit = Reap(
                prognameInit,
                self.zmat,
                initDisp.DispCart,
                self.options,
                initDisp.n_coord,
                eigsInit,
                indices,
                self.options.energyRegexInit,
                self.options.successRegexInit,
            )
            ReapObjInit.run()
            # os.chdir('..')
            print("GiraffeDirectory:")
            print(os.getcwd())
            # raise RuntimeError

            # nate
            p_en_arrayInit = ReapObjInit.p_en_array
            m_en_arrayInit = ReapObjInit.m_en_array
            ref_enInit = ReapObjInit.ref_en

            fcInit = ForceConstant(
                initDisp,
                p_en_arrayInit,
                m_en_arrayInit,
                ref_enInit,
                self.options,
                indices,
            )
            fcInit.run()
            print("Computed Force Constants:")
            print(fcInit.FC)

            # raise RuntimeError

        if os.path.exists(rootdir + "/fc.dat") or os.path.exists(rootdir + "/FCMFINAL"):
            f_read.run()
            f_conv = F_conv(
                f_read.FC_mat,
                s_vec,
                self.zmat,
                "internal",
                False,
                self.TED,
                self.options.units,
            )
            f_conv.run()
            F = f_conv.F
        else:
            F = fcInit.FC

        if self.options.coords != "ZMAT" and not initBool:
            F = np.dot(self.TED.Proj.T, np.dot(F, self.TED.Proj))
        if self.options.coords != "ZMAT":
            g_mat.G = np.dot(self.TED.Proj.T, np.dot(g_mat.G, self.TED.Proj))

        """
            Run the GF matrix method with the internal F-Matrix and computed G-Matrix!
        """
        print("Initial Frequencies:")
        init_GF = GF_Method(
            g_mat.G.copy(),
            F.copy(),
            self.options.tol,
            self.options.projTol,
            self.zmat,
            self.TED,
        )
        init_GF.run()

        """
            Now for the TED check.
        """
        self.G = np.dot(np.dot(LA.inv(init_GF.L), g_mat.G), LA.inv(init_GF.L).T)
        self.G[np.abs(self.G) < self.options.tol] = 0
        self.F = np.dot(np.dot(init_GF.L.T, F), init_GF.L)
        self.F[np.abs(self.F) < self.options.tol] = 0

        print("TED Frequencies:")
        TED_GF = GF_Method(
            self.G, self.F, self.options.tol, self.options.projTol, self.zmat, self.TED
        )
        TED_GF.run()

        # if os.path.exists(rootdir + '/Intder'):
        # shutil.rmtree(rootdir + '/Intder')
        S = TED_GF.S
        initial_fc = TED_GF.eig_v
        eigs = len(S)

        """ Coming soon, machinery that computes a diagnostic to aid in off-diagonal force constant
        computation priority!! """

        algo = Algorithm(eigs, S, initial_fc, self.options)
        algo.run()
        print(algo.indices)

        """
            Recompute the B-Tensors to match the final geometry, 
            then generate the displacements.
        """

        s_vec = s_vectors(self.zmat, self.options, self.zmat.variableDictionaryFinal)
        s_vec.run(self.zmat.CartesiansFinal, False)
        transdisp = TransDisp(
            s_vec,
            self.zmat,
            self.options.disp,
            init_GF.L,
            True,
            self.options.dispTol,
            self.TED,
            self.options,
            algo.indices,
            GF=TED_GF,
        )
        transdisp.run()
        # nate
        eigs = transdisp.eigs
        p_disp = transdisp.p_disp
        m_disp = transdisp.m_disp

        if self.options.dispCheck:
            raise RuntimeError

        """
            The displacements have been generated, now we have to run them!
        """

        prog = self.options.program
        progname = prog.split("@")[0]

        if self.options.calc:
            Dir_obj = DirectoryTree(
                progname,
                self.zmat,
                transdisp,
                self.options.cartInsert,
                p_disp,
                m_disp,
                self.options,
                algo.indices,
                "template.dat",
            )
            Dir_obj.run()
            os.chdir(rootdir + "/Disps")
            dispList = []
            for i in os.listdir(rootdir + "/Disps"):
                dispList.append(i)
            print("printing disp list")
            print(dispList)

            """
                This code generates the submit script for the displacements, 
                submits an array, then checks if all jobs have finished every 
                10 seconds. :D
            """

            v_template = vulcan_template(self.options, len(dispList), progname, prog)
            out = v_template.run()
            with open("displacements.sh", "w") as file:
                file.write(out)

            pipe = subprocess.PIPE

            process = subprocess.run(
                "qsub displacements.sh", stdout=pipe, stderr=pipe, shell=True
            )
            self.outRegex = re.compile(r"Your\s*job\-array\s*(\d*)")
            self.job_id = int(re.search(self.outRegex, str(process.stdout)).group(1))

            self.jobFinRegex = re.compile(r"taskid")
            while True:
                qacct_proc = subprocess.run(
                    ["qacct", "-j", str(self.job_id)], stdout=pipe, stderr=pipe
                )
                qacct_string = str(qacct_proc.stdout)
                job_match = re.findall(self.jobFinRegex, qacct_string)
                if len(job_match) == len(dispList):
                    break
                time.sleep(10)

            output = str(process.stdout)
            error = str(process.stderr)

        """
            After this point, all of the jobs will have finished, and its time 
            to reap the energies as well as checking for sucesses on all of 
            the jobs
        """
        if not self.options.calc:
            os.chdir("Disps")
        Reap_obj = Reap(
            progname,
            self.zmat,
            transdisp.DispCart,
            self.options,
            transdisp.n_coord,
            eigs,
            algo.indices,
            self.options.energyRegex,
            self.options.successRegex,
        )
        Reap_obj.run()
        # os.chdir('..')

        # nate
        p_en_array = Reap_obj.p_en_array
        m_en_array = Reap_obj.m_en_array
        ref_en = Reap_obj.ref_en

        """
            Compute the force constants here, currently can only do diagonal 
            force constants
        """
        # fc = ForceConstant(transdisp, Reap_obj.energiesDict)

        # nate
        fc = ForceConstant(
            transdisp, p_en_array, m_en_array, ref_en, self.options, algo.indices
        )
        fc.run()
        print("Computed Force Constants:")
        print(fc.FC)

        """
            I will need to make a method just for reading in force constants, 
            for now the diagonals will do.
        """
        # self.F = np.diag(fc.FC)
        # nate
        self.F = fc.FC
        """
            Recompute the G-matrix with the new geometry, and then transform 
            the G-matrix using the lower level of theory eigenvalue matrix. 
            This will not fully diagonalize the G-matrix if a different 
            geometry is used between the two.
        """
        g_mat = G_Matrix(self.zmat, s_vec, self.options)
        g_mat.run()
        if self.options.coords != "ZMAT":
            g_mat.G = np.dot(self.TED.Proj.T, np.dot(g_mat.G, self.TED.Proj))
        self.G = np.dot(np.dot(transdisp.eig_inv, g_mat.G), transdisp.eig_inv.T)
        self.G[np.abs(self.G) < self.options.tol] = 0

        """
            Final GF Matrix run
        """
        print("Final Frequencies:")
        Final_GF = GF_Method(
            self.G, self.F, self.options.tol, self.options.projTol, self.zmat, self.TED
        )
        Final_GF.run()

        """
            This code below is a rudimentary table of the TED for the final 
            frequencies. Actually right now it uses the initial L-matrix, 
            which may need to be modified by the final L-matrix.
        """
        print("////////////////////////////////////////////")
        print("//{:^40s}//".format(" Final TED"))
        print("////////////////////////////////////////////")
        self.TED.run(np.dot(init_GF.L, Final_GF.L), Final_GF.Freq)

        """
            This code prints out the frequencies in order of energy as well 
            as the ZPVE in several different units.
        """
        print(
            "Final ZPVE in: "
            + "{:6.2f}".format(np.sum(Final_GF.Freq) / 2)
            + " (cm^-1) "
            + "{:6.2f}".format(0.5 * np.sum(Final_GF.Freq) / 349.7550881133)
            + " (kcal mol^-1) "
            + "{:6.2f}".format(0.5 * np.sum(Final_GF.Freq) / 219474.6313708)
            + " (hartrees) "
        )

        """
            This code converts the force constants back into cartesian 
            coordinates and writes out an "output.default.hess" file, which 
            is of the same format as FCMFINAL of CFOUR.

            One further note. I will likely have to convert my force constants 
            within the F_Conv back to internal coordinates first.
        """
        self.F = np.dot(np.dot(transdisp.eig_inv.T, self.F), transdisp.eig_inv)
        if self.options.coords != "ZMAT":
            self.F = np.dot(self.TED.Proj, np.dot(self.F, self.TED.Proj.T))
        cart_conv = F_conv(
            self.F, s_vec, self.zmat, "cartesian", True, self.TED, self.options.units
        )
        cart_conv.run()

        t2 = time.time()

        # print(Final_GF.Freq - init_GF.Freq)
        print("This program took " + str(t2 - t1) + " seconds to run.")
