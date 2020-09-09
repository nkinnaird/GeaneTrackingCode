
localsetup=$5/localProducts_gm2_v8_04_00_prof/setup

echo "ls $localsetup : `ls -lh $localsetup`"
echo "ls $5 : `ls $5`"
echo "ls $5/local* : `ls $5/local*`"

[ ! -f $localsetup -o ! -r $localsetup ] && echo -e "\nUnable to access local setup file $localsetup\n" && exit

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup
setup ifdhc v1_6_2 -z /cvmfs/fermilab.opensciencegrid.org/products/common/db

echo "Here is the your environment in this job: " > job_output_${CLUSTER}.${PROCESS}.log
env >> job_output_${CLUSTER}.${PROCESS}.log

if [ "${PROCESS}" == "0" ]; then
    echo ${JOBSUBJOBID} > ${JOBSUBJOBID}
    ifdh cp -D ${JOBSUBJOBID} $3
fi

#Manipulate file names
startingFile=$4
startingFileRoot=$(echo $startingFile | cut -d. -f1)
startingFileNum=$(echo $startingFile | cut -d. -f2)
startingFileExt=$(echo $startingFile | cut -d. -f3)

#Write file to pass for ifdh cp command
rm -f fileCopyList.txt && touch fileCopyList.txt
echo $startingFileRoot.${PROCESS}.$startingFileExt . >> fileCopyList.txt


#Copy over input files
#Moving this further down causes proxyCert error?!
echo "fileCopyList.txt:"
cat fileCopyList.txt
echo "ifdh cp -f fileCopyList.txt"
ifdh cp -f fileCopyList.txt

for file in `ls | grep $(basename $startingFileRoot)`; do
	if [ ! -s "$file" ]; then
		echo "$file is empty. Deleting..."
		rm -f $file
	fi
done

source /cvmfs/gm2.opensciencegrid.org/prod7/g-2/setup

# export PRODUCTS="/gm2/app/users/nkinnaird/localProducts:${PRODUCTS}"

source $localsetup

#####IMPORTANT!!!#######
#first time and then every time you change something in your local development area, 
#you should install it before running from grid
#########################
#mrb i --generator ninja

#need to setup local products before running
. mrb slp

gm2 -c $1.fcl -n $2 -s $(basename $startingFileRoot).*

mv strawTrackLeftRightGEANE.log strawTrackLeftRightGEANE_${CLUSTER}.${PROCESS}.log
mv strawLeftRightGEANETracks.root strawLeftRightGEANETracks_${CLUSTER}.${PROCESS}.root
mv GEANEPlots.root GEANEPlots_${CLUSTER}.${PROCESS}.root

ifdh cp -D job_output_${CLUSTER}.${PROCESS}.log $3/logs
ifdh cp -D strawTrackLeftRightGEANE_${CLUSTER}.${PROCESS}.log $3/logs
ifdh cp -D strawLeftRightGEANETracks_${CLUSTER}.${PROCESS}.root $3/data
ifdh cp -D GEANEPlots_${CLUSTER}.${PROCESS}.root $3/data


