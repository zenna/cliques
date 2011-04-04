#!/bin/bash

for (( i=3; i<=5; i++))
do
for (( j=1; j<=5; j++))
do
  n1=$((5*$i))
  n2=$((4*$j))
  n3=3
  k1=$(($n1-1))
  k2=$(($i*($n2-1)))
  k3=$(($i*2*$j*($n3-1)))
  graph_file="${n1}_${n2}_${n3}_${k1}_${k2}_${k3}"
  ./random_graph --random.n1=${n1} --random.n2=${n2} --random.n3=${n3} --random.k1=${k1} --random.k2=${k2} --random.k3=${k3} ${graph_file}.bin > ${graph_file}.txt
  ./master.py ${graph_file}.bin $((n1*n2*n3))
done
done
