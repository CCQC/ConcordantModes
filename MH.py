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
from MixedHessian.Final_intder_input import intder_final
from MixedHessian.F_Convert          import F_conv
from MixedHessian.F_Read             import F_Read
from MixedHessian.ForceConstant      import ForceConstant
from MixedHessian.G_Matrix           import G_Matrix
from MixedHessian.GF_Method          import GF_Method
from MixedHessian.Gen_Final_Intder   import Final_Intder
from MixedHessian.GrabEig            import GrabEig
from MixedHessian.GrabSym            import GrabSym
from MixedHessian.init_intder_input  import intder_init
from MixedHessian.Reap               import Reap
from MixedHessian.s_vectors          import s_vectors
from MixedHessian.TED_intder_input   import intder_100
from MixedHessian.TED                import TED
from MixedHessian.TransDisp          import TransDisp
from MixedHessian.vulcan_template    import vulcan_template
from MixedHessian.ZMAT_parse         import ZMAT

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
        f_conv = F_conv(np.identity(len(f_read.FC_mat)), f_read.FC_mat, "cart", s_vec, self.zmat)
        f_conv.run()


        """
            Compute G-Matrix
        """
        L = np.identity(len(s_vec.B))
        g_mat = G_Matrix(L, self.zmat, s_vec)
        g_mat.run()
        # raise RuntimeError

        """
            Run the GF matrix method with the internal F-Matrix and computed G-Matrix!
        """
        self.tol = 1e-14
        print("Initial Frequencies:")
        # init_GF = GF_Method(np.identity(len(g_mat.G)),g_mat.G.copy(),f_conv.F.copy(),1e-14,"Initial")
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
        # TED_GF = GF_Method(init_GF.L,g_mat.G.copy(),f_conv.F.copy(),1e-14,"TED")
        TED_GF = GF_Method(self.G,self.F,self.tol)
        TED_GF.run()
        # print('Transformed Force Constants')
        # print(np.dot(np.dot(np.transpose(init_GF.L),f_conv.F),init_GF.L))

        """
            Next block, generate and then run intder initial
        """
        if os.path.exists(rootdir + '/Intder'):
            shutil.rmtree(rootdir + '/Intder')
        # os.mkdir('Intder')
        # shutil.copy('FCMFINAL','Intder/file15')
        # os.chdir('Intder')
        # self.intder_init = intder_init(self.zmat)
        # self.intder_init.run()
        # os.system('/home/vulcan/mel64643/bin/MixedHessian/Temporary_Scripts/INTDER < intder.inp > intder.out')
        # # raise RuntimeError
       
        
        # """
            # Some post processing of INTDER initial to generate the SALCS for intder 100
        # """
        # GrabE = GrabEig()
        # GrabE.run()
        # ted = TED()
        # ted.run()
        
        # """
            # move some of the INTDER initial files, then generate and run the 100 INTDER job
        # """
        # shutil.move('intder.inp','init_intder.inp')
        # shutil.move('intder.out','init_intder.out')
        # shutil.copy('file15','init_file15')
        # shutil.copy('file16','init_file16')
        
        # self.intder_100 = intder_100(self.zmat)
        # self.intder_100.run()
        # os.system('/home/vulcan/mel64643/bin/MixedHessian/Temporary_Scripts/INTDER < intder.inp > intder.out')
        # # raise RuntimeError
        
        # """
            # script that grabs the Normal internal coordinate values
        # """
        # grabSym = GrabSym()
        # grabSym.run()
        
        # """
            # Move some INTDER files around, then move onto the next thing
        # """
        # shutil.copy('intder.out','../intderTEDcheck.out')
        # shutil.copy('eigen.csv','../')
        # shutil.copy('symVariables.csv','../')
        # shutil.move('intder.inp','100_intder.inp')
        # shutil.move('intder.out','100_intder.out')
        # shutil.move('file15','100_file15')
        # shutil.copy('file16','100_file16')
        
        # os.chdir('..')
        
        """
            This should replace the first part of the mathematica intdif script
        """
        # transdisp = TransDisp(s_vec,self.zmat,self.options.rdisp,self.options.adisp,init_GF.L)
        transdisp = TransDisp(s_vec,self.zmat,self.options.rdisp,self.options.adisp,init_GF.L)
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
        Reap_obj = Reap(progname,self.zmat,transdisp.DispCart)
        Reap_obj.run()
        os.chdir('..')

        """
            This section will compute the force constants using a python script
        """
        # fc = ForceConstant(self.zmat ,self.options.rdisp ,self.options.adisp, Reap_obj.energiesDict)
        fc = ForceConstant(self.zmat, transdisp, Reap_obj.energiesDict)
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
        self.G = np.dot(np.dot(transdisp.eig_inv,g_mat.G),np.transpose(transdisp.eig_inv))
        self.G[np.abs(self.G) < self.tol] = 0

        """
            Final GF Matrix run
        """
        
        print("Final Frequencies:")
        Final_GF = GF_Method(self.G,self.F,self.tol)
        Final_GF.run()
        
        """
            And now we scoot back over to the Intder directory to run the final hessian!
        """
        # os.chdir('Intder')
        # FinalIntder = intder_final(self.zmat)
        # FinalIntder.run()
        # Fin_Int = Final_Intder(fc.FC)
        # Fin_Int.run()
        # shutil.copy('intder.out','../MixedHessOutput.dat')
        # os.chdir('..')
        
        t2 = time.time()
        
        print('This program took ' + str(t2-t1) + ' seconds to run.')
