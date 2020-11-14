import fileinput
import os
import re
import shutil


# Alright, it's time to generalize this for multiple program inputs
    
class DirectoryTree(object):
    def __init__(self,progname,zmat,disps,insertionIndex):
        self.progname = progname
        self.zmat = zmat
        self.insertionIndex = insertionIndex
        self.disps = disps # This should be the 'TransDisp' object

    def Make_Input(self,data,dispp,n_at,at,index):	
        """
            I'm going to generalize this so that the user just puts in the insertion index,
            and the geometry is placed there. Then the user can use any program/input file that this
            package is coded for.
        """
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
        if os.path.exists(root+'initden.dat'):
            init = True
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
        os.chdir('..')
        for i in range(len(self.disps.DispCart)-1):
            j = i % 2
            k = i // 2
            os.mkdir(str(i+2))
            os.chdir("./" + str(i+2))
            if j == 0:
                data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_plus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
            if j == 1:
                data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_minus'],str(n_atoms),self.zmat.atomList,self.insertionIndex)
            with open(inp,'w') as file:
               file.writelines(data)
            data = data_buff.copy()
            if init:
                shutil.copy('../../initden.dat','.')
            os.chdir('..')

