import fileinput
import os
import re
import shutil
import numpy as np


# Alright, it's time to generalize this for multiple program inputs
    
class DirectoryTree(object):
    def __init__(self,progname,zmat,disps,insertionIndex,p_disp,m_disp,options,indices):
    # def __init__(self,progname,zmat,disps,insertionIndex,dispSym):
        # self.dispSym = dispSym
        self.progname = progname
        self.zmat = zmat
        self.insertionIndex = insertionIndex
        self.disps = disps # This should be the 'TransDisp' object
        #nate
        self.p_disp = p_disp
        self.m_disp = m_disp
        self.options = options
        self.indices = indices 
    def Make_Input(self,data,dispp,n_at,at,index):	 
        space = ' '
        if self.progname == 'cfour':
            space = ''
        if index == -1:
            print('The user needs to specify a different value for the cartInsert keyword.')
            raise RuntimeError
        else:
            for i in range(int(n_at)):
                data.insert(index+i, space + at[i] + "{:16.10f}".format(dispp[i][0]) + "{:16.10f}".format(dispp[i][1]) + "{:16.10f}".format(dispp[i][2]) + '\n')
        
        return data

    def run(self):
        disp_dict = {}
        disp_dict['disp'] = []
        disp_dict['geom'] = []

        root = os.getcwd()
        packagepath = os.path.realpath(__file__)
        packagepath = packagepath[:-len('/DirectoryTree.py')]

        progname = self.progname

        n_atoms = len(self.zmat.atomList)
        
        if progname == 'molpro' or progname == 'psi4' or progname =='cfour':
            with open('template.dat','r') as file:
                data = file.readlines()
        else:
            print('Specified program not supported: ' + progname)
            raise RuntimeError
        
        init = False
        genbas = False
        if os.path.exists(root+'/initden.dat'):
            init = True
        if os.path.exists(root+'/GENBAS'):
            genbas = True
        data_buff = data.copy()
        if os.path.exists(os.getcwd() + '/Disps'):
            shutil.rmtree('Disps',ignore_errors=True)
        os.mkdir('Disps')
        os.chdir('./Disps')
        os.mkdir("1")
        os.chdir("./1")
        data = self.Make_Input(data,self.disps.DispCart['ref'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
        inp = ''
        if self.progname == 'cfour':
            inp = 'ZMAT'
        else:
            inp = 'input.dat'

        with open(inp, 'w') as file:
            file.writelines(data)
        data = data_buff.copy()
        if init:
            shutil.copy('../../initden.dat','.')
        if genbas:
            shutil.copy('../../GENBAS','.')
        os.chdir('..')

        """
            Not including the reference input, this function generates the directories for the displacement
            jobs and copies in the input file data. Following this, these jobs are ready to be submitted to the queue.  
        """ 
        
    #def create_directory(self,direc,data):
    #    os.mkdir(str(direc))
    #    os.chdir("./" + str(direc))          
    #    print(data)
    #    file = open(inp,'w')
    #    file.writelines(data)
    #    if init:
    #        shutil.copy('../../initden.dat','.')
    #    if genbas:
    #        shutil.copy('../../GENBAS','.')
    #    os.chdir('..')
    #    file.close()
        p_disp = self.p_disp
        m_disp = self.m_disp
        a = p_disp.shape[0]
        indices = self.indices 
        direc = 2 
       
        """
            This loop makes calls to the Make_Input function for storing the cartesian displacements as a variable,
            upon which it is written to a standard input file after the create_directory function creates a directory
            within the Disps directory.
        """
         
        #for index in indices: 
        #    i,j = index[0], index[1] 
        #    p_data = self.Make_Input(data,p_disp[i,j],str(n_atoms),self.zmat.atomList,self.insertionIndex)  

        #    self.create_directory(direc,p_data)
        #    m_data = self.Make_Input(data,m_disp[i,j],str(n_atoms),self.zmat.atomList,self.insertionIndex) 

        #    self.create_directory(direc+1,m_data)
        #    direc += 2              
        
        for index in indices: 
            i,j = index[0], index[1] 
            p_data = self.Make_Input(data,p_disp[i,j],str(n_atoms),self.zmat.atomList,self.insertionIndex)  

            #create_directory(direc,p_data)
            os.mkdir(str(direc))
            os.chdir("./" + str(direc))          
            with open(inp,'w') as file:
                file.writelines(p_data)
            data = data_buff.copy()
            if init:
                shutil.copy('../../initden.dat','.')
            if genbas:
                shutil.copy('../../GENBAS','.')
            os.chdir('..')

            m_data = self.Make_Input(data,m_disp[i,j],str(n_atoms),self.zmat.atomList,self.insertionIndex) 
            os.mkdir(str(direc+1))
            os.chdir("./" + str(direc+1))          
            with open(inp,'w') as file:
                file.writelines(m_data)
            data = data_buff.copy()
            if init:
                shutil.copy('../../initden.dat','.')
            if genbas:
                shutil.copy('../../GENBAS','.')
            os.chdir('..')

            #create_directory(direc+1,m_data)
            #data = data_buff.copy()
            direc += 2              

        #old diagonal-only code
        
        #""" I need to restructure this so that it iterates over only running jobs, not the duplicates by symmetry. """
        #Sum = 0
        ## self.dispSym = self.dispSym.astype(int)
        ## for i in range(len(self.disps.DispCart)-np.sum(self.dispSym)-1):
        #for i in range(len(self.disps.DispCart)-1):
        #    j = Sum % 2
        #    k = Sum // 2
        #    os.mkdir(str(i+2))
        #    os.chdir("./" + str(i+2))
        #    if j == 0:
        #        data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_plus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
        #    elif j == 1:
        #        data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_minus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
        #    # if self.dispSym[k] == 0:
        #        # if j == 0:
        #            # data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_plus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
        #        # if j == 1:
        #            # data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_minus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
        #    # elif self.dispSym[k] == 1:
        #        # data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_plus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
        #    # Sum += self.dispSym[k]
        #    Sum += 1
        #    with open(inp,'w') as file:
        #        file.writelines(data)
        #    data = data_buff.copy()
        #    if init:
        #        shutil.copy('../../initden.dat','.')
        #    if genbas:
        #        shutil.copy('../../GENBAS','.')
        #    os.chdir('..')
