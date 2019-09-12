#!/bin/bash



for f in {1..100}; do
  result="$(find "$f" -maxdepth 1 -type f -name '*.csv' -printf .)"
  #echo "Folder '${f}' contains ${#result} *.csv files."
  #echo ${#result}
  if [[ ${#result} -eq 0 ]];
  then  
    echo "${f}"
  fi
done




