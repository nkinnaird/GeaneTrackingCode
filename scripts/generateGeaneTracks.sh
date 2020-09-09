#!/bin/sh
# Usage: generateGeaneTracks.sh [FCL filename] [# events] [# jobs] [FILE TO READ-IN]
# Common Usage: ./generateGeaneTracks.sh RunGeane -1 1000 /pnfs/path/to/events/gasgunfile.0.root

[ ! -r $PWD/submit-geane-tracks.sh ] && echo "Unable to access submit-geane-tracks.sh at $PWD !" && exit

#three levels up
localdir=$(dirname $(dirname $(dirname $PWD)))

#create a place to put your output. This creates a directory in the /pnfs area
#which is based on your Fermilab user name, and the current date and time. 

echo "Add folder tag information to end of Tracking directory:"
read folderNameTag

NOW=$(date +"%F-%H-%M-%S")
SCRATCH_DIR=/pnfs/GM2/scratch/users/${USER}/Tracking/${NOW}-$folderNameTag

echo "Scratch directory: ${SCRATCH_DIR}"

mkdir ${SCRATCH_DIR}
mkdir ${SCRATCH_DIR}/logs
mkdir ${SCRATCH_DIR}/data

chmod -R g+w ${SCRATCH_DIR}

#enter the FHiCL file you want to use for the submission
MAINFCLNAME=${1:-ProductionMuPlusMuonGasGun}
echo "Main FCL: ${MAINFCLNAME}.fcl"

cp $localdir/localProducts_gm2_v8_04_00_prof/gm2tracker/v8_04_00/fcl/${MAINFCLNAME}.fcl ${SCRATCH_DIR}
cp $localdir/localProducts_gm2_v8_04_00_prof/gm2tracker/v8_04_00/fcl/geaneFitParams.fcl ${SCRATCH_DIR}
cp $localdir/localProducts_gm2_v8_04_00_prof/gm2tracker/v8_04_00/fcl/strawCommonParams.fcl ${SCRATCH_DIR}


#How many events per job and how many jobs
NEVTSPERJOB=${2:-1}
echo "# of events per job: $NEVTSPERJOB"
NJOBS=${3:-1}
echo "# of jobs: $NJOBS"
INPUTFIRSTFILE=$4
echo "Input file: $INPUTFIRSTFILE"
# INPUTNUMFILES=$5
# echo "Number of input files: $INPUTNUMFILES"

#This submits the job to the grid using local release:
# echo jobsub_submit -N ${NJOBS} -G gm2 --OS=SL6  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=1h --role=Analysis file://$PWD/submit-geane-tracks.sh ${MAINFCLNAME} ${NEVTSPERJOB} ${SCRATCH_DIR} ${INPUTFIRSTFILE} $localdir # ${INPUTNUMFILES} $localdir
jobsub_submit -N ${NJOBS} -G gm2 --OS=SL6  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=2h --role=Analysis file://$PWD/submit-geane-tracks.sh ${MAINFCLNAME} ${NEVTSPERJOB} ${SCRATCH_DIR} ${INPUTFIRSTFILE} $localdir # ${INPUTNUMFILES} $localdir

