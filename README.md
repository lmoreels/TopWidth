# TopWidth

## Set up CMSSW to get correct ROOT version, ...

export SCRAM_ARCH=slc6_amd64_gcc491

cmsrel CMSSW_7_4_2

cd CMSSW_7_4_2/src

cmsenv

git cms-init

## Get TopTreeProducer from git

git clone https://github.com/TopBrussels/TopTreeProducer TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_74X

cd src

make

cd ../../..

## Get TopTreeAnalysisBase from git

git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/

git checkout CMSSW_74X

make

cd ../../

##Get private code directory from git

Make a new local repository and check the box to create a readme file.

git clone https://github.com/lmoreels/TopWidth TopWidth

git checkout master
