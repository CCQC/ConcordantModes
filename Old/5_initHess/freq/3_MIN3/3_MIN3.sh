#!/bin/sh
#$ -q gen6.q
#$ -N MIN6_Mixed_H
#$ -S /bin/sh
#$ -cwd

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
echo "    Submitted using:   submit -N MIN6_Mixed_H -n 4 gen6.q cfour@2.0+mpi"
echo "***********************************************************************"


# Load the requested Cfour module file
vulcan load cfour@2.0+mpi

scratch=$TMPDIR/$USER/$JOB_ID

# Set variables
export OMP_NUM_THREADS=1
export NSLOTS=4
prefix=/opt/vulcan/opt/vulcan/linux-x86_64/intel-16.0.1/cfour-2.0-hjkzrfqmrmfvexnyvlji4zmpovcjpspy

# Copy job data
cp $SGE_O_WORKDIR/ZMAT $scratch/ZMAT
cp $prefix/basis/GENBAS $scratch
cp $prefix/basis/ECPDATA $scratch
if [ -e JAINDX ]; then cp JAINDX $scratch ; fi
if [ -e JOBARC ]; then cp JOBARC $scratch ; fi
if [ -e FCMINT ]; then cp FCMINT $scratch ; fi
if [ -e GENBAS ]; then cp GENBAS $scratch ; fi
if [ -e ECPDATA ]; then cp ECPDATA $scratch ; fi
if [ -e OPTARC ]; then cp OPTARC $scratch ; fi
if [ -e ISOTOPES ]; then cp ISOTOPES $scratch ; fi
if [ -e ISOMASS ]; then cp ISOMASS $scratch ; fi
if [ -e initden.dat ]; then cp initden.dat $scratch ; fi

echo " Running cfour on `hostname`"
echo " Running calculation..."

cd $scratch
xcfour >& $SGE_O_WORKDIR/output.dat
xja2fja
/opt/scripts/cfour2avogadro $SGE_O_WORKDIR/output.dat

echo " Saving data and cleaning up..."
if [ -e ZMATnew ]; then cp -f ZMATnew $SGE_O_WORKDIR/ZMATnew ; fi

# Create a job data archive file
tar --transform "s,^,Job_Data_$JOB_ID/," -vcf $SGE_O_WORKDIR/Job_Data_$JOB_ID.tar OPTARC FCMINT FCMFINAL ZMATnew JMOL.plot JOBARC JAINDX FJOBARC DIPDER HESSIAN MOLDEN NEWMOS AVOGADROplot.log den.dat
if [ -e zmat001 ]; then tar --transform "s,^,Job_Data_$JOB_ID/," -vrf $SGE_O_WORKDIR/Job_Data_$JOB_ID.tar zmat* ; fi
gzip $SGE_O_WORKDIR/Job_Data_$JOB_ID.tar

echo " Job complete on `hostname`."
