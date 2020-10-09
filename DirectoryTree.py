import fileinput
import os
import re
import shutil


# Alright, it's time to generalize this for multiple program inputs
    
class DirectoryTree(object):
    def __init__(self,progname,zmat,disps):
        self.progname = progname
        self.zmat = zmat
        self.disps = disps # This should be the 'TransDisp' object

    def Make_Input(self,data,dispp,n_at,at,prog):	
        """
            Modify this later so that f-strings are used for the geometry... but for now it works!
        """
        self.psi4GeomRegex = re.compile(r"molecule\s+\{")
        self.molproGeomRegex = re.compile(r"geometry\s+=\s+\{")
        self.endGeomRegex = re.compile(r"}")
        if prog == 'molpro':
            ind1 = 0
            for i in range(len(data)):
                if re.search(self.molproGeomRegex,data[i]):
                    ind1 = i
            for i in range(len(data)-ind1-1):
                if re.search(self.endGeomRegex,data[i+ind1+1]):
                    ind2 = i+ind1+1
                    break
            for i in range(int(n_at)):
                data.insert(ind2+i, ' ' + at[i] + "{:16.10f}".format(dispp[i][0]) + "{:16.10f}".format(dispp[i][1]) + "{:16.10f}".format(dispp[i][2]) + '\n')
        elif prog == 'psi4':
            ind1 = 0
            for i in range(len(data)):
                if re.search(self.psi4GeomRegex,data[i]):
                    ind1 = i
            for i in range(len(data)-ind1-1):
                if re.search(self.endGeomRegex,data[i+ind1+1]):
                    ind2 = i+ind1+1
                    break
            for i in range(int(n_at)):
                data.insert(ind2+i, '  ' + at[i] + "{:16.10f}".format(dispp[i][0]) + "{:16.10f}".format(dispp[i][1]) + "{:16.10f}".format(dispp[i][2]) + '\n')
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
        
        if progname == 'molpro' or progname == 'psi4':
            with open('template.dat','r') as file:
                data = file.readlines()
        else:
            print('Specified program not supported: ' + progname)
            raise RuntimeError
        
        data_buff = data.copy()
        if os.path.exists(os.getcwd() + '/Disps'):
            shutil.rmtree('Disps',ignore_errors=True)
        os.mkdir('Disps')
        os.chdir('./Disps')
        os.mkdir("1")
        os.chdir("./1")
        data = self.Make_Input(data,self.disps.DispCart['ref'],str(n_atoms),self.zmat.atomList,progname)
        with open('input.dat','w') as file:
            file.writelines(data)
        data = data_buff.copy()
        os.chdir('..')
        for i in range(len(self.disps.DispCart)-1):
            j = i % 2
            k = i // 2
            os.mkdir(str(i+2))
            os.chdir("./" + str(i+2))
            if j == 0:
                data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_plus'],str(n_atoms),self.zmat.atomList,progname)
            if j == 1:
                data = self.Make_Input(data,self.disps.DispCart[str(k+1)+'_minus'],str(n_atoms),self.zmat.atomList,progname)
            with open('input.dat','w') as file:
               file.writelines(data)
            data = data_buff.copy()
            os.chdir('..')

