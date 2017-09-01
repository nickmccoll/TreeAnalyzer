#!/bin/bash

cmd=$1
outputdir=$2
CMSSWVERS=$3

workdir=`pwd`

echo `hostname`
echo "${_CONDOR_SCRATCH_DIR}"
echo "workdir: $workdir"
echo "args: $*"
tar -xf data.tar.gz --warning=no-timestamp
rm data.tar.gz
ls -l -a

source /cvmfs/cms.cern.ch/cmsset_default.sh
SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 project CMSSW ${CMSSWVERS}`
cd CMSSW_8_0_27/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
echo "CMSSW: "$CMSSW_BASE

cd $workdir
eval $cmd

### Now that the cmsRun is over, there is one or more root files created
echo "List all root files = "
ls *.root
echo "List all files"
ls 
echo "*******************************************"
for FILE in *.root
do
    if [[ "$outputdir" =~ ^/eos/uscms/.* ]]; then
        copypath=`echo ${outputdir} | sed "s:/eos/uscms::"`
        xrdcp -np ${FILE} root://cmseos:1094/${copypath}/
        rm ${FILE}
    fi
done

exit $status
