import shutil
import os
rootdir = os.getcwd()
# print(rootdir)
# print(os.path.realpath(__file__))

packagepath = os.path.realpath(__file__)
packagepath = packagepath[:-6]

# os.chdir(packagepath + '/subprocess_Scripts')
os.chdir(packagepath)
print('Package path:')
print(os.getcwd())
os.chdir(rootdir)
print('Working path:')
print(os.getcwd())

if os.path.exists(rootdir + '/subprocess_Scripts'):
    shutil.rmtree(rootdir + '/subprocess_Scripts')
if os.path.exists(rootdir + '/templates'):
    shutil.rmtree(rootdir + '/templates')



shutil.copytree(packagepath + '/subprocess_Scripts',rootdir + '/subprocess_Scripts')
shutil.copytree(packagepath + '/templates',rootdir + '/templates')

# Parse the output to get all pertinent ZMAT info
if os.path.exists(rootdir + '/zmatFiles'):
    shutil.rmtree(rootdir + '/zmatFiles')
os.system('python subprocess_Scripts/ZMAT_interp.py')

# Next block, generate and then run intder initial
if os.path.exists(rootdir + '/Intder_init'):
    shutil.rmtree(rootdir + '/Intder_init')
os.mkdir('Intder_init')
shutil.copy('templates/intder_init_template.dat','Intder_init')
shutil.copy('FCMFINAL','Intder_init/file15')
os.chdir('Intder_init')
os.system('python ../subprocess_Scripts/init_intder_input.py')
os.system(packagepath + '/subprocess_Scripts/INTDER < intder.inp > intder.out')

# Some post processing of INTDER initial to generate the SALCS for intder 100
os.system('python ../subprocess_Scripts/GrabEig.py')
os.system('python ../subprocess_Scripts/TED.py')

# delete some of the INTDER initial files and bring in the 100 files
os.remove('intder.inp')
os.remove('intder.out')

shutil.copy('../templates/intder_100_template.dat','.')
os.system('python ../subprocess_Scripts/100_intder_input.py')
os.system(packagepath + '/subprocess_Scripts/INTDER < intder.inp > intder.out')
shutil.copy('intder.out','../intderTEDcheck.out')
shutil.copy('eigen.csv','../')
os.remove('intder.inp')
os.remove('intder.out')
os.remove('file15')

os.chdir('..')


# Now for the Mathematica portion of our adventure

if os.path.exists(rootdir + '/mma'):
    shutil.rmtree(rootdir + '/mma')
os.mkdir('mma')
os.chdir('mma')
shutil.copy('../subprocess_Scripts/Intdif2008.m','.')
shutil.copy('../eigen.csv','.')
os.system('python ../subprocess_Scripts/GenLoad.py')
os.system('./Load.wls')




# if os.path.exists(rootdir + '/zmatFiles'):
    # shutil.rmtree(rootdir + '/zmatFiles')
# if os.path.exists(rootdir + '/subprocess_Scripts'):
    # shutil.rmtree(rootdir + '/subprocess_Scripts')
# if os.path.exists(rootdir + '/templates'):
    # shutil.rmtree(rootdir + '/templates')
