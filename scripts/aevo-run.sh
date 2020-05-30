#!/bin/bash

final=0

#----------------------------------------------------------------------#

if [[ -f state.json ]]
then
  gen=$(grep generation state.json | sed 's/^.*: \([0-9]*\),/\1/g')
else
  ./find-groundstate.py input.json > ga.out
  gen=0
fi

for i in $(seq $gen $final)
do

  echo -n "Running generation ${i} ... "

  cd generation$(printf "%05d" $i)
  ../eval-fitness.sh
  cd ..

  ./find-groundstate.py input.json --state state.json >> ga.out

  echo "done."

done
