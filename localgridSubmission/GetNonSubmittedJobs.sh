#!/bin/bash

path="./"

if [ $# -eq 1 ]
then path=$1
fi

echo "Checking if there are submit scripts without output log (*.sh.o*) ..."

for f in "$path"*.sh; do
  x=`ls -d $path$f* | wc -l`;
  if [ $x == 1 ]
  then echo $f
  fi
done
