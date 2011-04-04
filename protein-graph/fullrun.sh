#!/bin/bash
# fullrun.sh
# Invokes gromacs.sh & graph_convert.sh and deals with directory structure
# First argument must be *.pdb
# Yun William Yu (c) 2010

# Find script path
script_path="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
bin_path=`dirname "$script_path"`
echo This is the location of the scripts: $bin_path

name=${1%.pdb}		# Strip extension off pdb file

# Temporary working directory for all the loose files
origdir=`pwd`
tmpdir=${name}.tmp
mkdir $tmpdir
cp ${name}.pdb $tmpdir

cd $tmpdir

#Rename the mdp files in question
for file in ${bin_path}/*.mdp; do
  #echo $file
  sed "s/\[insert_name_here\]/${name}/" $file > ${name}`basename $file`
done

${bin_path}/gromacs-emonly.sh ${name}.pdb
if [ "$?" -ne "0" ] ; then
	echo gromacs.sh exited with error code $?
	exit 1
fi

cat ${name}_MD.pdb | grep -v SOL | grep -v NA+ | grep -v CL- > ${name}_MD_mod.pdb
${bin_path}/graph_convert.sh  ${name}_MD_mod.pdb

cd ${name}_MD_mod.graphtmp
mkdir ../${name}_matlab
cp *.out saltbridges.txt ../${name}_matlab
cd ..
cp ${name}allH.top ${name}_matlab

cp -a ${bin_path}/pdb2adj.orig/* ${name}_matlab

cd ${name}_matlab
echo cd `pwd` >> matlab_run.m
echo "adjacency('cov.out','${name}allH.top','covtable.txt','hbonds.out','saltbridges.txt','hphobes5A.out','hphobes8A.out');" >> matlab_run.m
echo exit >> matlab_run.m
matlab -nodesktop -nodisplay -nosplash  < matlab_run.m > matlab_run.out

: <<'END'
END
