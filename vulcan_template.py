import numpy as np
from numpy.linalg import inv
from numpy import linalg as LA


"""
    This file just stores the vulcan qsub script
"""

class vulcan_template(object):
    def __init__(self,options,job_num,progname,prog):
        self.progdict = {
            "molpro": "molpro -n $NSLOTS --nouse-logfile --no-xml-output -o output.dat input.dat",
            "psi4": "psi4 -n $NSLOTS"
        }
        self.odict = {                               
            'q':        options.queue,                      
            'nslots':   options.nslots,    
            'jarray':   '1-{}'.format(job_num), 
            'progname': progname,               
            'prog':     prog,                   
            'tc':       str(job_num),           
            'cline':    self.progdict[progname]      
        }
# This can be inserted back in if the sync keyword is sorted
#$ -sync y
        self.vulcan_template = """#!/bin/sh
#$ -q {q}
#$ -N MixedHess
#$ -S /bin/sh
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

{cline}"""
    def run(self):
        return self.vulcan_template.format(**self.odict)
