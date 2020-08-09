import fileinput
import os
import shutil

def Make_Input(data,dispp,n_at,at):	
    data.insert(7,' ' + n_at + " \n" )
    for i in range(int(n_at)):
        data.insert(9+i , ' ' + at[i] + dispp[i])
    return data
    
if __name__ == "__main__":
    disp_dict = {}
    disp_dict['disp'] = []
    disp_dict['geom'] = []

    root = os.getcwd()
    packagepath = os.path.realpath(__file__)
    packagepath = packagepath[:-len('/DirectoryTree.py')]

    os.chdir('zmatFiles')

    with open("atomList",'r') as file:
        atoms = file.read()
    atoms = atoms.split('\n')
    
    os.chdir(root)

    n_atoms = len(atoms)
    dispCart = "dispcart" 
    
    os.chdir('mma')
    
    with open(dispCart,'r') as file:
        disp = file.readlines()
   
    os.chdir(root)
    
    with open("basis",'r') as file:
        b = file.readlines()
    
    # print(b[0][:-1])
    basis = b[0]

    # disp = disp[:-1]
    n_disp = int((len(disp))/(n_atoms+1))
    
    for i in range(n_disp):
        disp_dict['disp'].append(i + 1)
        disp_dict['geom'].append([])
        for j in range(n_atoms):
            disp_dict['geom'][i].append(disp[i*(n_atoms + 1) + j + 1])
    
    os.chdir(packagepath + '/../templates')

    with open('input.dat','r') as file:
        data = file.readlines()
    
    data[10] = 'basis=' + basis
    print(data[10])
    os.chdir(root)

    data_buff = data.copy()
    if os.path.exists(os.getcwd() + '/Disps'):
        shutil.rmtree('Disps')
    os.mkdir('Disps')
    os.chdir('./Disps')
    for i in range(n_disp):
        os.mkdir(str(i+1))
        os.chdir("./" + str(i+1))
        data = Make_Input(data,disp_dict['geom'][i],str(n_atoms),atoms)
        with open('input.dat','w') as file:
           file.writelines(data)
        data = data_buff.copy()
        os.chdir('..')

