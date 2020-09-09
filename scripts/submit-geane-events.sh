
localsetup=$4/localProducts_gm2_v8_04_00_prof/setup

echo "ls $localsetup : `ls -lh $localsetup`"
echo "ls $4 : `ls $4`"
echo "ls $4/local* : `ls $4/local*`"

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

cp $4/localProducts_gm2_v8_04_00_prof/gm2ringsim/v8_04_00/fcl/$1.fcl .
echo "source.firstRun: 10${CLUSTER}"  > runinfo
echo "source.firstSubRun: ${PROCESS}" >> runinfo
cat $1.fcl runinfo > theFCL.fcl

gm2 -c theFCL.fcl -n $2 


 mv gm2ringsim_mdc0-geane.log gm2ringsim_mdc0-geane_${CLUSTER}.${PROCESS}.log
 mv gm2ringsim_muon_gas_gun-geane.root gm2ringsim_muon_gas_gun-geane_${CLUSTER}.${PROCESS}.root

 ifdh cp -D job_output_${CLUSTER}.${PROCESS}.log $3/logs
 ifdh cp -D gm2ringsim_mdc0-geane_${CLUSTER}.${PROCESS}.log $3/logs
 ifdh cp -D gm2ringsim_muon_gas_gun-geane_${CLUSTER}.${PROCESS}.root $3/data

# mv gm2ringsim_mdc2-geane.log gm2ringsim_mdc2-geane_${CLUSTER}.${PROCESS}.log
# mv gm2ringsim_mdc2-geane.root gm2ringsim_mdc2-geane_${CLUSTER}.${PROCESS}.root

# ifdh cp -D job_output_${CLUSTER}.${PROCESS}.log $3/logs
# ifdh cp -D gm2ringsim_mdc2-geane_${CLUSTER}.${PROCESS}.log $3/logs
# ifdh cp -D gm2ringsim_mdc2-geane_${CLUSTER}.${PROCESS}.root $3/data
