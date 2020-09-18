import fileinput
import os
import shutil


# Alright, it's time to generalize this for multiple program inputs
    
class DirectoryTree(object):
    def __init__(self,progname,basis,charge,spin,zmat,disps):
        self.progname = progname
        self.basis = basis
        self.charge = charge
        self.spin = spin
        self.zmat = zmat
        self.disps = disps # This should be the 'TransDisp' object

    def Make_Input(self,data,dispp,n_at,at,prog):	
        """
            Modify this later so that f-strings are used for the geometry... but for now it works!
        """
        if prog == 'molpro':
            data.insert(7,' ' + n_at + " \n" )
            for i in range(int(n_at)):
                data.insert(9+i , ' ' + at[i] + "{:16.10f}".format(dispp[i][0]) + "{:16.10f}".format(dispp[i][1]) + "{:16.10f}".format(dispp[i][2]) + '\n')
                # data.insert(9+i , ' ' + at[i] + ' ' + str(dispp[i][0]) + ' ' + str(dispp[i][1]) + ' ' + str(dispp[i][2]) + '\n')
        elif prog == 'psi4':
            data[5] = '  ' + str(self.charge) + ' ' + str(self.spin) + '\n'
            if self.spin != 1:
                data[11] = '  reference uhf\n'
            for i in range(int(n_at)):
                data.insert(6+i , '  ' + at[i] + "{:16.10f}".format(dispp[i][0]) + "{:16.10f}".format(dispp[i][1]) + "{:16.10f}".format(dispp[i][2]) + '\n')
                # data.insert(6+i , '  ' + at[i] + ' ' + str(dispp[i][0]) + ' ' + str(dispp[i][1]) + ' ' + str(dispp[i][2]) + '\n')
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
        # dispCart = "dispcart" 
        
        # os.chdir('mma')
        
        # with open(dispCart,'r') as file:
            # disp = file.readlines()
   
        # os.chdir(root)
        
        basis = self.basis

        # n_disp = int((len(disp))/(n_atoms+1))
        
        # for i in range(n_disp):
            # disp_dict['disp'].append(i + 1)
            # disp_dict['geom'].append([])
            # for j in range(n_atoms):
                # disp_dict['geom'][i].append(disp[i*(n_atoms + 1) + j + 1])
        
        os.chdir(packagepath + '/templates')

        basis_string = ''
        if progname == 'molpro':
            with open('input.dat','r') as file:
                data = file.readlines()
        elif progname == 'psi4':
            with open('psi4_template.dat','r') as file:
                data = file.readlines()
            basis_string += '  '
        else:
            print('Specified program not supported: ' + progname)
            raise RuntimeError
        
        # raise RuntimeError
        basis_string += 'basis=' + basis + '\n'
        data[10] = basis_string
        # print(data[10])
        os.chdir(root)

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

