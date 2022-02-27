#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export HOME=/afs/cern.ch/user/a/athachay
export X509_USER_PROXY=/afs/cern.ch/user/a/athachay/private/.proxy/x509up_u134523
cd /afs/cern.ch/work/a/athachay/private/l1egamma/isolut/CMSSW_12_0_2/src/L1EGCalibrations/Isolation/test/Condor/JobsForSteps/v2
eval `scramv1 runtime -sh`
set -x
TMPDIR=`mktemp -d`
cd $TMPDIR
cp  /afs/cern.ch/work/a/athachay/private/l1egamma/isolut/CMSSW_12_0_2/src/L1EGCalibrations/Isolation/test/applyIsolation.exe .
./applyIsolation.exe /afs/cern.ch/work/a/athachay/private/l1egamma/isolut/CMSSW_12_0_2/src/L1EGCalibrations/Isolation/test/Condor/JobsForSteps/v2/ParList_ApplyIsolaion.dat
if [ $? -eq 0 ]; then 
    mv * /afs/cern.ch/work/a/athachay/private/l1egamma/isolut/CMSSW_12_0_2/src/L1EGCalibrations/Isolation/test/Condor/JobsForSteps/v2
else
    echo FAIL
fi
