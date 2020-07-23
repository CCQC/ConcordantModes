#!/bin/sh
#$ -q gen4.q
#$ -N WaterTest
#$ -S /bin/sh
#$ -cwd
#$ -t 1-7
#$ -tc 4

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
echo "    Submitted using:   submit -N WaterTest --array 1-7 gen4.q molpro@2010.1.67+mpi"
echo "***********************************************************************"

# cd into individual task directory
cd $SGE_O_WORKDIR/$SGE_TASK_ID
vulcan load molpro@2010.1.67+mpi

export NSLOTS=4

molpro -n $NSLOTS --nouse-logfile --no-xml-output -o output.dat input.dat
