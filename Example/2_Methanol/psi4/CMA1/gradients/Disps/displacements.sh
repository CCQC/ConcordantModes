#!/bin/sh
#$ -q gen5.q,gen4.q,gen6.q
#$ -N Concordant
#$ -S /bin/sh
#$ -cwd
#$ -t 1-25
#$ -tc 25

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
echo "    Submitted using:   Concordant Modes Approach "
echo "***********************************************************************"

# cd into individual task directory
cd $SGE_O_WORKDIR/$SGE_TASK_ID
vulcan load psi4@master

export NSLOTS=1

psi4 -n $NSLOTS