# TopWidth

## Set up CMSSW to get correct ROOT version, ...

export SCRAM_ARCH=slc6_amd64_gcc491

cmsrel CMSSW_7_4_2

cd CMSSW_7_4_2/src

cmsenv

git cms-init

## Get TopTreeProducer from git

git clone https://github.com/TopBrussels/TopTreeProducer TopTreeProducer

cd TopTreeProducer/

git checkout CMSSW_74X

cd src

make

cd ../..

scram b -j8
(I did this after I got all the packages, but then there are problems with TopTreeAnalysisBase.)

## Get TopTreeAnalysisBase from git

git clone https://github.com/TopBrussels/TopTreeAnalysisBase TopTreeAnalysisBase/

cd TopTreeAnalysisBase/

git checkout CMSSW_74X

make

cd ..

## Get private code directory from git

--Make a new local repository on the website and check the box to create a readme file.
Now you can clone this repository in different ways:
If you use the https link, you will need to enter your username and password each time you push, pull or fetch.
If you use the ssh link, you need to enter the passphrase of your ssh key.--

git clone ssh://git@github.com/lmoreels/TopWidth TopWidth

git checkout master


## Compile

source compile.sh

## Run

./testAnalyser
