#!/bin/bash

#PBS -q localgrid
#PBS -l walltime=08:00:00

#source /user/lmoreels/.bashrc
#source $VO_CMS_SW_DIR/cmsset_default.sh
# setting up your code and your env
cd /user/lmoreels/CMSSW_7_6_5/src/ 
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd TopBrussels/TopWidth
