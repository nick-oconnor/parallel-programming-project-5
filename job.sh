#!/bin/bash
rm results.txt
for j in 64 16 4 1
do
  echo "$j tasks with $((64/$j)) threads on 64 nodes"
  srun --partition=small --time=15 --overcommit --runjob-opts="--mapping TEDCBA" --nodes=64 --ntasks-per-node=$j ./a.out $((64/$j)) >> results.txt
  if [ $? -ne 0 ]
    then exit 0
  fi
  echo "DONE!"
done
for i in 128 256
do
  for j in 64 16 4 1
  do
    echo "$j tasks with $((64/$j)) threads on $i nodes"
    srun --partition=medium --time=15 --overcommit --runjob-opts="--mapping TEDCBA" --nodes=$i --ntasks-per-node=$j ./a.out $((64/$j)) >>results.txt
    if [ $? -ne 0 ]
      then exit 0
    fi
    echo "DONE!"
  done
done
