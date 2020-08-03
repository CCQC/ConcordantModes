import time
import sys
import shutil
import os
import subprocess

t1 = time.time()

rootdir = os.getcwd()

packagepath = os.path.realpath(__file__)
packagepath = packagepath[:-len('/MH.py')]

# First things first. Read in options from the MHopt.py file



# Parse the output to get all pertinent ZMAT info
if os.path.exists(rootdir + '/zmatFiles'):
    shutil.rmtree(rootdir + '/zmatFiles')
os.system('python ' + packagepath + '/subprocess_Scripts/ZMAT_interp.py')

# Next block, generate and then run intder initial
if os.path.exists(rootdir + '/Intder'):
    shutil.rmtree(rootdir + '/Intder')
os.mkdir('Intder')
shutil.copy('FCMFINAL','Intder/file15')
os.chdir('Intder')
os.system('python ' + packagepath + '/subprocess_Scripts/init_intder_input.py')
os.system('/home/vulcan/mel64643/bin/MixedHessian/subprocess_Scripts/INTDER < intder.inp > intder.out')

# Some post processing of INTDER initial to generate the SALCS for intder 100
os.system('python ' + packagepath + '/subprocess_Scripts/GrabEig.py')
os.system('python ' + packagepath + '/subprocess_Scripts/TED.py')

# move some of the INTDER initial files, then generate and run the 100 INTDER job
shutil.move('intder.inp','init_intder.inp')
shutil.move('intder.out','init_intder.out')

os.system('python ' + packagepath + '/subprocess_Scripts/100_intder_input.py')
os.system('/home/vulcan/mel64643/bin/MixedHessian/subprocess_Scripts/INTDER < intder.inp > intder.out')

# script that grabs the symmetry internal coordinate values
os.system('python ' + packagepath + '/subprocess_Scripts/GrabSym.py')


# Move some INTDER files around, then move onto the next thing
shutil.copy('intder.out','../intderTEDcheck.out')
shutil.copy('eigen.csv','../')
shutil.copy('symVariables.csv','../')
shutil.move('intder.inp','100_intder.inp')
shutil.move('intder.out','100_intder.out')
shutil.move('file15','old_file15')

os.chdir('..')


# Now for the Mathematica portion of our adventure
if os.path.exists(rootdir + '/mma'):
    shutil.rmtree(rootdir + '/mma')
os.mkdir('mma')
os.chdir('mma')
# Unfortunately, Intdif must be in the same directory as the mathematica script...for now
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

q = 'gen4.q,gen5.q,gen6.q'
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

# Remove scratch files in cwd()
# os.remove('eigen.csv')
# os.remove('symVariables.csv')

# if os.path.exists(rootdir + '/zmatFiles'):
    # shutil.rmtree(rootdir + '/zmatFiles')
# if os.path.exists(rootdir + '/mma'):
    # shutil.rmtree(rootdir + '/mma')
# if os.path.exists(rootdir + '/Intder'):
    # shutil.rmtree(rootdir + '/Intder')

t2 = time.time()

print('This program took ' + str(t2-t1) + ' seconds to run.')
