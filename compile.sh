#!/bin/bash

# if there is one arg compile only this .cc
if [[ -n $1 ]]
then
    if [[ $1 == *".cc" ]]
    then ccfile=$1
    else ccfile=$1".cc"
    fi

    ofile=`echo $ccfile |sed 's/\.cc$//g'`
    echo "compiling : " $ccfile ", executible name: " $ofile
    if [[ "$OSTYPE" == "linux-gnu" ]];
    then g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent76 -l TopTreeAna76 -l TopWidthAna76 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
    elif [[ "$OSTYPE" == "darwin"* ]];
    then g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent76 -l TopTreeAna76 -l TopWidthAna76 -l MLP -l TreePlayer -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
    fi


# if there is no arg compile all .cc
else
    for ccfile in ./*.cc
    do
        if [[ $ccfile == *"Dict"* ]] ||  [[ $ccfile == *"_local"* ]]
        then continue
        fi
        
        ofile=`echo $ccfile |sed 's/\.cc$//g'`
        echo "compiling : " $ccfile ", executible name: " $ofile
        if [[ "$OSTYPE" == "linux-gnu" ]]
        then g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent76 -l TopTreeAna76 -l TopWidthAna76 -l MLP -l TreePlayer -l TMVA -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
        elif [[ "$OSTYPE" == "darwin"* ]]
        then g++ -g -std=c++11 -L ~/lib -L . -L .. -I ./ -I ../ -l TopTreeAnaContent76 -l TopTreeAna76 -l TopWidthAna76 -l MLP -l TreePlayer -l XMLIO -I `root-config --incdir` `root-config --libs` $ccfile -o $ofile
        fi
    done
fi
