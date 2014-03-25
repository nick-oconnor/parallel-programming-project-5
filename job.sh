#!/bin/bash

for i in 128 256
  do
  for j in 1 2 4
    do
    srun --partition=medium --time=10 --overcommit --runjob-opts="--mapping TEDCBA" --nodes=$i --ntasks=$(($i*$j)) ./MPI3
    echo $i $j
  done
done
