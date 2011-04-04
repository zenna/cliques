#!/bin/bash
# This file uses first to convert a pdb file into a graph
# Copyright (c) 2010 Yun William Yu

name=${1%.pdb}	# Strip extension off pdb file

# Temporary working directory for all the loose files
origdir=`pwd`
tmpdir=${name}.graphtmp
mkdir $tmpdir
cp ${name}.pdb $tmpdir
cd $tmpdir

# Remove solvent, if necessary
cat ${name}.pdb | grep -v SOL | grep -v NA+  | grep -v NA+ | sed -e 's/O1  GLY/O   GLY/' -e 's/O2  GLY/OXT GLY/' -e 's/CD  ILE/CD1 ILE/' -e 's/O1/O /g' -e 's/O2 /OXT/g' > ${name}-noSOL.pdb
g_editconf -f ${name}-noSOL.pdb -o ${name}-noSOL2.pdb
mv ${name}-noSOL2.pdb ${name}-noSOL.pdb
mv ${name}-noSOL.pdb ${name}.pdb

FIRST_ROOT=/home/ywy09/Downloads/FIRST-6.2.1-bin-32-gcc3.4.3-O3 /home/ywy09/Downloads/FIRST-6.2.1-bin-32-gcc3.4.3-O3/FIRST ./${name}.pdb -c 2.0 -covout -E 0.01 -H 3 -hbout -o 3 -phout -v 3  -non

# Extract salt bridges
cat ${name}_results.txt | grep SALT | sed "s/COMPUTING ENERGY FOR SALT BRIDGE (D-H - A):[^0-9]*[0-9]*//" > ./saltbridges.txt
cat ./saltbridges.txt

# FIRST doesn't understand metal ions (yet), so this following couple of lines add in support for MG++ ions
cat ${name}_results.txt | sed -e '/| MG|MG|/d' -e '/MG/!d' | awk '{ print $1, $6 , 5 }' >> cov.out


# Create yet another directory to hold the files with different energy cutoff
mkdir 8A
cp ${name}.pdb 8A
cd 8A
FIRST_ROOT=/home/ywy09/Downloads/FIRST-6.2.1-bin-32-gcc3.4.3-O3 /home/ywy09/Downloads/FIRST-6.2.1-bin-32-gcc3.4.3-O3/FIRST ./${name}.pdb -c 4.5 -covout -E -0.01 -H 3 -hbout -o 3 -phout -v 3  -non

# Rename appropriate hphobes files
cp hphobes.out ../hphobes8A.out
cd ..
cp hphobes.out hphobes5A.out
cd ..

