#!/bin/bash 

if [[ -n $1 ]] #check if variable is not empty
then
    if [[ $1 == "test" ]]
    then
        cd test
        for f in ./submit*.sh
        do
            qsub -q express $f
        done
        cd -
    elif [[ $1 == "sendnow" ]]
    then
        cd output
        for f in ./submit*.sh
        do
            qsub $f
        done
        cd -
    fi

else
    cd output
    outputFile="listSubmit.txt"
    for f in ./submit*.sh
    do
        echo "qsub -q localgrid $f" >> $outputFile
    done
    big-submission $outputFile
    cd -

fi
