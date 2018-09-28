#!/usr/bin/env bash

arr=($*)
echo "mode ${arr[0]}"


MODE="${arr[0]}" 

unset arr[0]

# join rest of the array

ARGUMENTS=""

for var in "${arr[@]}"
do
  ARGUMENTS="$ARGUMENTS ${var}"
  # do something on $var
done

echo $ARGUMENTS


Rscript -e "source(file = system.file(\"exec/clustmut_$MODE.R\", package = \"clustMut\"))" $ARGUMENTS



