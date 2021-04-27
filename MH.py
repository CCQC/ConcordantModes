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
from MixedHessian.DirectoryTree      import DirectoryTree
from MixedHessian.F_Convert          import F_conv
from MixedHessian.F_Read             import F_Read
from MixedHessian.ForceConstant      import ForceConstant
from MixedHessian.G_Matrix           import G_Matrix
from MixedHessian.GF_Method          import GF_Method
from MixedHessian.Reap               import Reap
from MixedHessian.s_vectors          import s_vectors
from MixedHessian.TED                import TED
from MixedHessian.TransDisp          import TransDisp
from MixedHessian.vulcan_template    import vulcan_template
from MixedHessian.ZMAT               import ZMAT
from MixedHessian.int2cart           import int2cart

class MixedHessian(object):
    def __init__(self, options):
        self.options = options
        """ 
            This constant is from:
            https://physics.nist.gov/cgi-bin/cuu/Value?hr 
            If this link dies, find the new link on NIST
            for the Hartree to Joule conversion and pop
            it in there.
            There is a standard uncertainty of 0.0000000000085
            to this value.
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
        s_vec = s_vectors(self.zmat,self.options)
        s_vec.run(self.zmat.CartesiansInit,True)

        self.TED = TED(s_vec.Proj,self.zmat)
        
        """
            Compute G-Matrix
        """
        g_mat = G_Matrix(self.zmat, s_vec, self.options)
        g_mat.run()
        
        G = g_mat.G.copy()

        """
            Read in FC matrix in cartesians, then convert to internals.
        """
        if os.path.exists(rootdir + '/fc.dat'):
            f_read = F_Read("fc.dat")
        elif os.path.exists(rootdir + '/FCMFINAL'):
            f_read = F_Read("FCMFINAL")
        else:
            print("Need to specify the force constants!")
            raise RuntimeError
        f_read.run()
        f_conv = F_conv(f_read.FC_mat, s_vec, self.zmat, "internal", False, self.TED)
        f_conv.run()
        if self.options.coords != "ZMAT":
            f_conv.F = np.dot(self.TED.Proj.T,np.dot(f_conv.F,self.TED.Proj))
        if self.options.coords != "ZMAT":
            g_mat.G = np.dot(self.TED.Proj.T,np.dot(g_mat.G,self.TED.Proj))
        
        """
            Run the GF matrix method with the internal F-Matrix and computed G-Matrix!
        """
        print("Initial Frequencies:")
        init_GF = GF_Method(g_mat.G.copy(),f_conv.F.copy(),self.options.tol,self.options.projTol,self.zmat,self.TED)
        init_GF.run()

        """
            Now for the TED check.
        """
        self.G = np.dot(np.dot(LA.inv(init_GF.L),g_mat.G),LA.inv(init_GF.L).T)
        self.G[np.abs(self.G) < self.options.tol] = 0
        self.F = np.dot(np.dot(init_GF.L.T,f_conv.F),init_GF.L)
        self.F[np.abs(self.F) < self.options.tol] = 0
        
        print("TED Frequencies:")
        TED_GF = GF_Method(self.G,self.F,self.options.tol,self.options.projTol,self.zmat,self.TED)
        TED_GF.run()

        if os.path.exists(rootdir + '/Intder'):
            shutil.rmtree(rootdir + '/Intder')
        
        """
            Recompute the B-Tensors to match the final geometry, then generate the displacements.
        """
        s_vec = s_vectors(self.zmat,self.options)
        s_vec.run(self.zmat.CartesiansFinal,False)
        transdisp = TransDisp(s_vec,self.zmat,self.options.disp,init_GF.L,True,self.options.dispTol,self.TED,self.options)
        transdisp.run()
        
        
        if self.options.dispCheck:
            raise RuntimeError

        """
            The displacements have been generated, now we have to run them!
        """
        
        prog = self.options.program
        progname = prog.split('@')[0]
        
        if self.options.calc:
            Dir_obj = DirectoryTree(progname, self.zmat, transdisp, self.options.cartInsert, transdisp.dispSym)
            Dir_obj.run()
            os.chdir(rootdir + '/Disps')
            dispList = []
            for i in os.listdir(rootdir + '/Disps'):
                dispList.append(i)
            
            """
                This code generates the submit script for the displacements, submits an array, then checks if all jobs have finished every 10 seconds. :D
            """

            v_template = vulcan_template(self.options,len(dispList),progname,prog)
            out = v_template.run()
            with open('displacements.sh','w') as file:
                file.write(out)
            
            pipe = subprocess.PIPE
            
            process = subprocess.run('qsub displacements.sh', stdout=pipe, stderr=pipe, shell=True)
            self.outRegex = re.compile(r'Your\s*job\-array\s*(\d*)')
            self.job_id = int(re.search(self.outRegex,str(process.stdout)).group(1))
            self.jobFinRegex = re.compile(r'taskid')
            while(True):
                qacct_proc = subprocess.run(['qacct','-j',str(self.job_id)], stdout=pipe, stderr=pipe)
                qacct_string = str(qacct_proc.stdout)
                job_match = re.findall(self.jobFinRegex,qacct_string)
                if len(job_match) == len(dispList):
                    break
                time.sleep(10)

            output = str(process.stdout)
            error = str(process.stderr)  
        

        """
            After this point, all of the jobs will have finished, and its time to reap the energies
            as well as checking for sucesses on all of the jobs
        """
        if not self.options.calc:
            os.chdir("Disps")
        Reap_obj = Reap(progname,self.zmat,transdisp.DispCart,self.options,transdisp.n_coord)
        Reap_obj.run()
        os.chdir('..')

        """
            Compute the force constants here, currently can only do diagonal force constants
        """
        fc = ForceConstant(transdisp, Reap_obj.energiesDict)
        fc.run()
        print('Computed Force Constants:')
        print(fc.FC)

        """
            I will need to make a method just for reading in force constants, for now the diagonals will do.
        """
        self.F = np.diag(fc.FC)

        """
            Recompute the G-matrix with the new geometry, and then transform the G-matrix using the 
            lower level of theory eigenvalue matrix. This will not fully diagonalize the G-matrix
            if a different geometry is used between the two.
        """
        g_mat = G_Matrix(self.zmat, s_vec, self.options)
        g_mat.run()
        if self.options.coords != "ZMAT":
            g_mat.G = np.dot(self.TED.Proj.T,np.dot(g_mat.G,self.TED.Proj))
        self.G = np.dot(np.dot(transdisp.eig_inv,g_mat.G),transdisp.eig_inv.T)
        self.G[np.abs(self.G) < self.options.tol] = 0

        """
            Final GF Matrix run
        """
        print("Final Frequencies:")
        Final_GF = GF_Method(self.G,self.F,self.options.tol,self.options.projTol,self.zmat,self.TED)
        Final_GF.run()
        

        """
            This code prints out the frequencies in order of energy as well as the ZPVE in several different units.
        """
        print("Final ZPVE in: " + "{:6.2f}".format(np.sum(Final_GF.Freq)/2) + " (cm^-1) " + "{:6.2f}".format(0.5*np.sum(Final_GF.Freq)/349.7550881133) + " (kcal mol^-1) " \
                + "{:6.2f}".format(0.5*np.sum(Final_GF.Freq)/219474.6313708) + " (hartrees) ")
        
        """
            This code below is a rudimentary table of the TED for the final frequencies.
            Actually right now it uses the initial L-matrix, which may need to be modified
            by the final L-matrix.
        """
        self.TED.run(init_GF.L,Final_GF.Freq)
        
        """
            This code converts the force constants back into cartesian coordinates and writes out
            an "output.default.hess" file, which is of the same format as FCMFINAL of CFOUR.

            One further note. I will likely have to convert my force constants within the F_Conv back to internal coordinates first.
        """
        # cart_conv = F_conv(self.F, s_vec, self.zmat, "cartesian", True)
        # cart_conv.run()

        t2 = time.time()
        
        # print(Final_GF.Freq - init_GF.Freq)
        print('This program took ' + str(t2-t1) + ' seconds to run.')
