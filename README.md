# TopWidth

## Set up CMSSW

--To get correct versions of ROOT, python, ... (You need ROOT 6)--

export SCRAM_ARCH=slc6_amd64_gcc491

cmsrel CMSSW_7_4_2

cd CMSSW_7_4_2/src

cmsenv

git cms-init

## Get TopTreeProducer from git

--Make sure to add the 'TopBrussels' directory. Otherwise the compilation later on will fail.--

git clone https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer

cd TopBrussels/TopTreeProducer/

git checkout CMSSW_74X

cd src

make

cd ../../..

## Get TopTreeAnalysisBase from git

git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopBrussels/TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/

git checkout CMSSW_74X

make

cd ../..

## Compile CMSSW and TopBrussels packages

scram b -j16


## Get private code directory from git

--Make a new local repository on the website and check the box to create a readme file.
Now you can clone this repository in different ways:
If you use the https link, you will need to enter your username and password each time you push, pull or fetch.
If you use the ssh link, you need to enter the passphrase of your ssh key.--

git clone ssh://git@github.com/lmoreels/TopWidth TopWidth

git checkout master


## Compile and make executables of all macros

source compile.sh

## Run

./testAnalyser
