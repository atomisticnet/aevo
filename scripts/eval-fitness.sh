#!/bin/bash

for f in POSCAR-*
do
  id=$(echo $f | sed 's/POSCAR-\(.*\)\.vasp/\1/g')
  if ! [[ -f FITNESS.$id ]]
  then
      mg energy $f '{"O":-2, "Li":1, "Sb":5, "Ni":3.5}'
      mv $f.energy FITNESS.$id
  fi
done
