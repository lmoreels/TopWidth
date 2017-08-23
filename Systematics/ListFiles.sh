#!/bin/bash

mainDir="/user/lmoreels/CMSSW_8_0_27/src/TopBrussels/TopWidth/";
outputFile="list_syst.txt"
nFiles=0

#check if path is declared
if [ $# -eq 0 ]
then echo "ERROR: Expecting argument. Path not declared. See help." >&2;
exit 1;
fi

if [[ $1 == "--help" ]] || [[ $1 == "-h" ]]
then
  echo "Use cases for ListFiles.sh";
  echo "       ./ListFiles.sh [path]";
  echo "       ./ListFiles.sh [path] [output file name]";
  exit 1;
fi

path=$1

if [ $# -eq 2 ]
then outputFile=$2
fi

# Init output file
echo -n "" > $outputFile

for entry in "$mainDir$path/result_"*; do 
  ((nFiles++));
  echo "$entry" >> $outputFile;
done

echo "Found $nFiles files"
