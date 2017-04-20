#!/bin/bash

# CHECK IF THESE ARE CORRECT AND ADAPT WHEN NECESSARY
CMSSWversion="CMSSW_80X_v8"
globalTag="TTP-CMSSW_80X_v8--GT-80X_mcRun2_asymptotic_2016_TrancheIV_v8"


#check if input file is declared
if [ $# -eq 0 ]
then echo "ERROR: Expecting argument. Input file not declared." >&2;
exit 1;
fi

inputFile=$1

#check if input file exists
if [ ! -e $inputFile ]
then echo "ERROR: Input file not found." >&2;
exit 1;
else echo "Input file is: $inputFile"
fi

if [ $# -eq 3 ]
then CMSSWversion=$2; globalTag=$3
fi

#strip extension off filename to make output file
inputFileName=${inputFile%.*}
outputFile="Output_"$inputFileName".txt"
pnfsPath="/pnfs/iihe/cms/store/user/fblekman/TopTree/"$CMSSWversion"/"$globalTag"/"


#put new empty line at end of file if not already there
sed -i $inputFile -e '$a\'
#number of samples in input file
z=($(wc $inputFile))
nofSamples=${z[0]}
currentSample=0

#initialise output file
echo "" > $outputFile

#take sample names from input file and fetch directory
while read c; do

((currentSample++))
#echo -ne "Looking for sample $currentSample/$nofSamples : $c"\\r;
echo "Looking for sample $currentSample/$nofSamples : $c";
find $pnfsPath$c -maxdepth 3 -mindepth 3 -type d -exec echo {} \; >> $outputFile;
echo "" >> $outputFile;

if [ $currentSample == $nofSamples ]
then echo "Processed $nofSamples samples";
fi

done < $inputFile

