# TopWidth

## Set up CMSSW

To get correct versions of ROOT, python, ... (You need ROOT 6)

~~~
export SCRAM_ARCH=slc6_amd64_gcc493

cmsrel CMSSW_7_6_5

cd CMSSW_7_6_5/src

cmsenv
~~~

## Get TopTreeProducer from git

Make sure to add the 'TopBrussels' directory. Otherwise the compilation later on will fail.

~~~
git clone -b CMSSW_76X https://github.com/TopBrussels/TopTreeProducer TopBrussels/TopTreeProducer/

cd TopBrussels/TopTreeProducer/src/

make

cd ../../..
~~~

## Get TopTreeAnalysisBase from git

~~~
git clone -b CMSSW_76X https://github.com/TopBrussels/TopTreeAnalysisBase TopBrussels/TopTreeAnalysisBase/

cd TopBrussels/TopTreeAnalysisBase/

make

cd ../..
~~~

## Compile CMSSW and TopBrussels packages

This is necessary if you want to produce TopTrees. The 'make' in 'TopTreeProducer' only compiles the TRoot* objects. For TopTree production you need the analysers as well.

~~~
scram b -j 16
~~~

## Get private code directory from git

Make a new local repository on the website and check the box to create a readme file.
Now you can clone this repository in different ways:
If you use the https link, you will need to enter your username and password each time you push, pull or fetch.
If you use the ssh link, you need to enter the passphrase of your ssh key.

~~~
git clone ssh://git@github.com/lmoreels/TopWidth TopBrussels/TopWidth

cd TopBrussels/TopWidth

make
~~~


## Get the TopTreeProducer version that was used to make your TopTrees

Versions of 'TopTreeProducer' are 'tagged' when TopTrees are made, generally indicated with a '_vX' at the end of the branch name. If you do not use the corresponding tag, you can get nasty error messages/seg faults when running over the TopTrees.

~~~
cd ../../TopBrussels/TopTreeProducer/

git checkout -b CMSSW_76X_v4

make

cd ../..
~~~

Note that you need to return to the non-tagged version to make changes. Since this branch already exists, you just need to do the following in the 'TopTreeProducer' directory.

~~~
git checkout CMSSW_76X
~~~


## Compile and make executables of all macros

~~~
source compile.sh
~~~

## Compile one macro

~~~
source compile.sh NtupleAnalyser.cc
~~~

## Run

~~~
./NtupleAnalyser
~~~
