#!/bin/bash

# Final iteration:
final=5

#-----------------------------------------------------------------------
#                     Do not edit contents below
#-----------------------------------------------------------------------

if [[ -f state.json ]]
then
  gen=$(grep generation state.json | sed 's/^.*: \([0-9]*\),/\1/g')
else
  aevolution.py input.json > aevo.out 2> aevo.err || exit 1
  gen=0
fi

for i in $(seq $gen $final)
do
  echo -n "Running generation ${i} ... "
  cd generation$(printf "%05d" $i)
  ../eval-fitness.py
  cd ..
  aevolution.py input.json --state state.json >> aevo.out 2> aevo.err || exit 1
  echo "done."
done
