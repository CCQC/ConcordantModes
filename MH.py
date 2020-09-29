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
from MixedHessian.TransDisp          import TransDisp
from MixedHessian.vulcan_template    import vulcan_template
from MixedHessian.ZMAT_parse         import ZMAT
from MixedHessian.ZMAT_interp        import ZMAT_interp

class MixedHessian(object):
    def __init__(self, options):
        self.options = options

    def run(self):
        t1 = time.time()
        
        rootdir = os.getcwd()
        
        packagepath = os.path.realpath(__file__)
        packagepath = packagepath[:-len('/MH.py')]
        """
            Parse the output to get all pertinent ZMAT info
        """
        if os.path.exists(rootdir + '/zmatFiles'):
            shutil.rmtree(rootdir + '/zmatFiles')
        self.zmat = ZMAT()
        self.zmat.run()
        self.finalZmat = ZMAT_interp()
        self.finalZmat.run()
        # raise RuntimeError

        """
            This is a new addition, the hope is to compute s-vectors right off the bat!
        """
        s_vec = s_vectors(self.zmat)
        s_vec.run()

        """
            Read in FC matrix in cartesians, then convert to internals.
        """

        f_read = F_Read("FCMFINAL")
        f_read.run()
        f_conv = F_conv(f_read.FC_mat, s_vec, self.zmat)
        f_conv.run()


        """
            Compute G-Matrix
        """
        g_mat = G_Matrix(self.zmat, s_vec)
        g_mat.run()
        # raise RuntimeError

        """
            Run the GF matrix method with the internal F-Matrix and computed G-Matrix!
        """
        self.tol = 1e-14
        print("Initial Frequencies:")
        init_GF = GF_Method(g_mat.G.copy(),f_conv.F.copy(),self.tol)
        init_GF.run()
        # raise RuntimeError

        """
            Now for the TED check.
        """
        self.G = np.dot(np.dot(inv(init_GF.L),g_mat.G),np.transpose(inv(init_GF.L)))
        self.G[np.abs(self.G) < self.tol] = 0
        self.F = np.dot(np.dot(np.transpose(init_GF.L),f_conv.F),init_GF.L)
        self.F[np.abs(self.F) < self.tol] = 0
        
        print("TED Frequencies:")
        TED_GF = GF_Method(self.G,self.F,self.tol)
        TED_GF.run()

        if os.path.exists(rootdir + '/Intder'):
            shutil.rmtree(rootdir + '/Intder')
        
        """
            This should replace the first part of the mathematica intdif script
        """
        s_vec = s_vectors(self.finalZmat)
        s_vec.run()
        transdisp = TransDisp(s_vec,self.finalZmat,self.options.rdisp,self.options.adisp,init_GF.L)
        transdisp.run()
        

        """
            The displacements have been generated, now we have to run them!
        """
        
        prog = self.options.program
        progname = prog.split('@')[0]
        
        Dir_obj = DirectoryTree(progname, self.options.basis, self.options.charge, self.options.spin, self.zmat, transdisp)
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
        Reap_obj = Reap(progname,self.finalZmat,transdisp.DispCart)
        Reap_obj.run()
        os.chdir('..')

        """
            This section will compute the force constants using a python script
        """
        fc = ForceConstant(self.finalZmat, transdisp, Reap_obj.energiesDict)
        fc.run()

        """
            This is temporary to fit with the INTDER units. Delete when you make your own GF matrix method.
        """
        fc.FC *= 4.35974394
        print('Computed Force Constants:')
        print(fc.FC)

        """
            I will need to make a method just for reading in force constants, for now the diagonals will do.
        """
        
        self.F = np.diag(fc.FC)
        g_mat = G_Matrix(self.finalZmat, s_vec)
        g_mat.run()
        self.G = np.dot(np.dot(transdisp.eig_inv,g_mat.G),np.transpose(transdisp.eig_inv))
        self.G[np.abs(self.G) < self.tol] = 0

        """
            Final GF Matrix run
        """
        
        print("Final Frequencies:")
        Final_GF = GF_Method(self.G,self.F,self.tol)
        Final_GF.run()
        
        t2 = time.time()
        
        print('This program took ' + str(t2-t1) + ' seconds to run.')
