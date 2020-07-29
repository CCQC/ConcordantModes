import sys
import shutil
import os
import subprocess

rootdir = os.getcwd()
# print(rootdir)
# print(os.path.realpath(__file__))

packagepath = os.path.realpath(__file__)
packagepath = packagepath[:-len('/MH.py')]

# os.chdir(packagepath + '/subprocess_Scripts')
os.chdir(packagepath)
# print('Package path:')
# print(os.getcwd())
os.chdir(rootdir)
# print('Working path:')
# print(os.getcwd())

if os.path.exists(rootdir + '/subprocess_Scripts'):
    shutil.rmtree(rootdir + '/subprocess_Scripts')
if os.path.exists(rootdir + '/templates'):
    shutil.rmtree(rootdir + '/templates')



# shutil.copytree(packagepath + '/subprocess_Scripts',rootdir + '/subprocess_Scripts')
# shutil.copytree(packagepath + '/templates',rootdir + '/templates')

# Parse the output to get all pertinent ZMAT info
if os.path.exists(rootdir + '/zmatFiles'):
    shutil.rmtree(rootdir + '/zmatFiles')
os.system('python ' + packagepath + '/subprocess_Scripts/ZMAT_interp.py')

# Next block, generate and then run intder initial
if os.path.exists(rootdir + '/Intder'):
    shutil.rmtree(rootdir + '/Intder')
os.mkdir('Intder')
# shutil.copy(packagepath +  '/templates/intder_init_template.dat','.')
shutil.copy('FCMFINAL','Intder/file15')
os.chdir('Intder')
os.system('python ' + packagepath + '/subprocess_Scripts/init_intder_input.py')
os.system(packagepath + '/subprocess_Scripts/INTDER < intder.inp > intder.out')

# Some post processing of INTDER initial to generate the SALCS for intder 100
os.system('python ' + packagepath + '/subprocess_Scripts/GrabEig.py')
os.system('python ' + packagepath + '/subprocess_Scripts/TED.py')

# delete some of the INTDER initial files and bring in the 100 files
os.remove('intder.inp')
os.remove('intder.out')
# shutil.move('intder.inp','old_intder.inp')
# shutil.move('intder.out','old_intder.out')

# shutil.copy(packagepath + '/templates/intder_100_template.dat','../')
os.system('python ' + packagepath + '/subprocess_Scripts/100_intder_input.py')
os.system(packagepath + '/subprocess_Scripts/INTDER < intder.inp > intder.out')

# Insert script here that grabs the symmetry internal coordinate values
os.system('python ' + packagepath + '/subprocess_Scripts/GrabSym.py')


shutil.copy('intder.out','../intderTEDcheck.out')
shutil.copy('eigen.csv','../')
shutil.copy('symVariables.csv','../')
os.remove('intder.inp')
os.remove('intder.out')
os.remove('file15')

os.chdir('..')


# Now for the Mathematica portion of our adventure

if os.path.exists(rootdir + '/mma'):
    shutil.rmtree(rootdir + '/mma')
os.mkdir('mma')
os.chdir('mma')
shutil.copy(packagepath + '/subprocess_Scripts/Intdif2008.m','.')
shutil.copy('../eigen.csv','.')
shutil.copy('../symVariables.csv','.')
os.system('python ' + packagepath + '/subprocess_Scripts/GenLoad.py')
os.system('./Load.wls')
os.chdir('..')

# The displacements have been generated, now we have to run them!
os.system('python ' + packagepath + '/subprocess_Scripts/DirectoryTree.py')
os.chdir('Disps')
dispList = []
for i in os.listdir(rootdir + '/Disps'):
    dispList.append(i)

q = 'debug.q,gen4.q,gen5.q,gen6.q'
prog = 'molpro@2010.1.67+mpi'
progname = 'molpro'
job_num = len(dispList)

# This portion is highly cumbersome and will have to split out into another script eventually,
# but for now I just want my code to work. :D

vulcan_template = """#!/bin/sh
#$ -q {q}
#$ -N MixedHess
#$ -S /bin/sh
#$ -sync y
#$ -cwd
#$ -t {jarray}
#$ -tc {tc}

. /etc/profile.d/modules.sh

# Disable production of core dump files
ulimit -c 0

echo ""
echo "***********************************************************************"
echo " Starting job:"
echo ""
echo "    Name:              "$JOB_NAME
echo "    ID:                "$JOB_ID
echo "    Hostname:          "$HOSTNAME
echo "    Working directory: "$SGE_O_WORKDIR
echo ""
echo "    Submitted using:   MixedHessian "
echo "***********************************************************************"

# cd into individual task directory
cd $SGE_O_WORKDIR/$SGE_TASK_ID
vulcan load {prog}

export NSLOTS={nslots}

{cline}
"""

progdict = {
    "molpro": "molpro -n $NSLOTS --nouse-logfile --no-xml-output -o output.dat input.dat",
    "psi4": "psi4 -n $NSLOTS"
}

odict = {
    'q':        q,
    'nslots':   4,
    'jarray':   '1-{}'.format(job_num),
    'progname': progname,
    'prog':     prog,
    'tc':       str(job_num),
    'cline':    progdict[progname]
}
out = vulcan_template.format(**odict)

with open('displacements.sh','w') as file:
    file.write(out)
subprocess.call('qsub displacements.sh', stderr=sys.stdout.buffer, shell=True)
# print(os.getcwd())


# After this point, all of the jobs will have finished, and its time to reap the energies
# as well as checking for sucesses on all of the jobs
os.system('python ' + packagepath + '/subprocess_Scripts/Reap.py')

# Now to generate the e.m file, then to generate the force constants!
os.chdir('../mma')
os.system('python ' + packagepath + '/subprocess_Scripts/GenEM.py')
os.system('python ' + packagepath + '/subprocess_Scripts/GenFC.py')
os.system('./FC.wls')

# And now we scoot back over to the Intder directory to run the final hessian!
os.chdir('../Intder')
os.system('python ' + packagepath + '/subprocess_Scripts/Final_intder_input.py')
os.system('python ' + packagepath + '/subprocess_Scripts/Gen_Final_Intder.py')
shutil.copy('intder.out','../MixedHessOutput.dat')
os.chdir('..')

# Remove templates
# shutil.remove('intder_init_template.dat')
# shutil.remove('intder_100_template.dat')
os.remove('eigen.csv')
os.remove('symVariables.csv')

if os.path.exists(rootdir + '/zmatFiles'):
    shutil.rmtree(rootdir + '/zmatFiles')
if os.path.exists(rootdir + '/mma'):
    shutil.rmtree(rootdir + '/mma')
if os.path.exists(rootdir + '/Intder'):
    shutil.rmtree(rootdir + '/Intder')
