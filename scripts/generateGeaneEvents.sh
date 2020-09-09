#!/bin/sh
# Usage: generateGeaneEvents.sh [FCL filename] [# events] [# jobs]
# Common Usage: ./generateGeaneEvents.sh mdc0-geane 5000 1000

[ ! -r $PWD/submit-geane-events.sh ] && echo "Unable to access submit-geane-events.sh at $PWD !" && exit

#three levels up
localdir=$(dirname $(dirname $(dirname $PWD)))

#create a place to put your output. This creates a directory in the /pnfs area
#which is based on your Fermilab user name, and the current date and time. 

echo "Add folder tag information to end of Events directory:"
read folderNameTag

NOW=$(date +"%F-%H-%M-%S")
SCRATCH_DIR=/pnfs/GM2/scratch/users/${USER}/Events/${NOW}-$folderNameTag

echo "Scratch directory: ${SCRATCH_DIR}"

mkdir ${SCRATCH_DIR}
mkdir ${SCRATCH_DIR}/logs
mkdir ${SCRATCH_DIR}/data

chmod -R g+w ${SCRATCH_DIR}

#enter the FHiCL file you want to use for the submission
MAINFCLNAME=${1:-ProductionMuPlusMuonGasGun}
echo "Main FCL: ${MAINFCLNAME}.fcl"

cp $localdir/localProducts_gm2_v8_04_00_prof/gm2ringsim/v8_04_00/fcl/${MAINFCLNAME}.fcl ${SCRATCH_DIR}

#How many events per job and how many jobs
NEVTSPERJOB=${2:-1}
echo "# of events per job: $NEVTSPERJOB"
NJOBS=${3:-1}
echo "# of jobs: $NJOBS"

#This submits the job to the grid using local release:
jobsub_submit -N ${NJOBS} -G gm2 --OS=SL6  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --expected-lifetime=2h --role=Analysis file://$PWD/submit-geane-events.sh ${MAINFCLNAME} ${NEVTSPERJOB} ${SCRATCH_DIR} $localdir

