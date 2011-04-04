#!/bin/bash
# First argument must be *.pdb
# Copyright (c) 2009 Yun William Yu

name=${1%.pdb}		# Strip extension off pdb file

# hydrogens addition
# -f is the input pdb file, -o the output gromacs file, -p the topology file (more complete description of the molecule), 
# -i include file for topolty (output, posre.itp), -q is the output pdb file (clean.pdb), -n is the index file (clean.idx)
# -ff is the force field chosen (oplsaa is an all-atom force field)
pdb2gmx -f ${name}.pdb -o ${name}allH.gro -p ${name}allH.top -i -q -n -ff oplsaa >& 01output.pdb2gmx

if [ "$?" -ne "0" ] ; then
	cat 01output.pdb2gmx
	echo pdb2gmx exited with nonzero error code
	exit 1
fi

# box
# The molecule -f is placed in a box, the result is in the file -o ..., -bt determines thte box type, -d specifies the
# distance between the molecule and the box.
editconf -f ${name}allH.gro -o ${name}allHbox.gro -bt cubic -d 0.6 >& 02output.editconf_box

if [ "$?" -ne "0" ] ; then
	cat 02output.editconf_box
	echo editconf exited with nonzero error code
	exit 1
fi

# solvation
# -cp solute, -cs solvent (Simple Point Charge water (SPC) by default), -p topology file (input and output), -o structure output file 
genbox -cp ${name}allHbox.gro -p ${name}allH.top -o ${name}allHsol.gro -cs >& 03output.genbox

if [ "$?" -ne "0" ] ; then
	cat 03output.genbox
	echo genbox exited with nonzero error code
	exit 1
fi

# name change
mv ${name}allH.top ${name}allHsol.top &> 04output.mv
mv "#${name}allH.top.1#" ${name}allH.top &> 05output.mv

# NEUTRALIZATION
# creation of dummy mdp for preprocessing before genion
touch dummy.mdp >& 06output.touch
# preprocessing for genion
# The gromacs preprocessor reads a molecular topology file, checks the validity of the file, expands the
# topology from a molecular description to an atomic description.
# -f (input) grompp file with MD parameters; -po (output) grompp file with MD parameters; -c (input) structure file; -pp (output) topology file; -o (output) Run file; 
grompp -f dummy.mdp -po dummyout.mdp -c ${name}allHsol.gro -p ${name}allHsol.top -pp ${name}allHsol_pre.top -o ${name}allHsol.tpr >& 07output.grompp_genion

if [ "$?" -ne "0" ] ; then
	cat 07output.grompp_genion
	echo grompp exited with nonzero error code
	exit 1
fi

# Determines the number of positive and negative charges that need to be added in order to neutralize the system
#nn=0 ; np=0; 
#fl_charge=`cat 07output.grompp_genion | grep 'System has non-zero total charge'| sed 's/System has non-zero total charge: //' | sed 's/e+/*10^/' | bc`
#int_charge=$(awk -v var="$fl_charge" 'BEGIN { rounded = sprintf("%.0f", var); print rounded }')
#echo $int_charge
#if [ $int_charge -gt "0" ]; then nn=$int_charge ; fi
#if [ $int_charge -lt "0" ]; then np=$((-1*$int_charge)) ; fi
#echo Neg=$nn, Pos=$np
# ions addition
# This replaces the solvent by monoatomic ions
# -s run input file; -o (output) structure file; -p (input/output) topology file; -np number of positive ions; -pname name of the positive ion; -g log file; -pot pdb file; -nname name of the negative ions
#echo 12 | genion -s ${name}allHsol.tpr -o ${name}allHion.gro -p ${name}allHsol_pre.top -np $np -pname Na+ -g -pot ${name}allHion.pdb -nn $nn -nname Cl- >& 08output.genion
#echo 12 | genion -s ${name}allHsol.tpr -o ${name}allHion.gro -p ${name}allHsol_pre.top -np $np -pname Na+ -g ion.log -nn $nn -nname Cl- >& 08output.genion
#echo 12 | genion -s ${name}allHsol.tpr -o ${name}allHion.gro -p ${name}allHsol_pre.top  -neutral -pname Na+ -nname Cl- -g genion.log >& 08output.genion

# Consolidate the solvent molecules in the topology file. Otherwise genion adds too many ions.
filename=${name}allHsol_pre.top
sol_mol=0
cp $filename $filename-new
while read line; do
	sed "/$line/d" $filename-new > $filename-new2
	mv $filename-new2 $filename-new
	add=$(echo $line | sed 's/SOL//') 
	sol_mol=$(($sol_mol+$add))
done < <( tail -n $(( $(cat $filename | wc -l) - $(grep " molecules " $filename  -n | sed 's/:\[ molecules \]//') + 1 )) $filename | grep SOL )
echo "SOL          $sol_mol" >> $filename-new
mv $filename-new $filename


echo SOL | genion -s ${name}allHsol.tpr -o ${name}allHion.gro -p ${name}allHsol_pre.top -pname NA+ -g genion.log -pot ${name}allHion.pdb -nname CL- -conc 0.0001  -neutral >& 08output.genion

if [ "$?" -ne "0" ] ; then
	cat 08output.genion
	echo genion exited with nonzero error code
	exit 1
fi

if [ -f ${name}allHion.gro ] ; then
  echo ${name}allHion.gro exists already
else
  cp ${name}allHsol.gro ${name}allHion.gro
fi

# name change
mv ${name}allHsol_pre.top ${name}allHion.top >& 09output.mv3
mv "#${name}allHsol_pre.top.1#" ${name}allHsol_pre.top >& 10output.mv4


# ENERGY MINIMIZATION
# preprocessing for EM with steepest descent
# The gromacs preprocessor reads a molecular topology file, checks the validity of the file, expands the
# topology from a molecular description to an atomic description.
# -f (input) grompp file with MD parameters; -po (output) grompp file with MD parameters; -c (input) structure file; -pp (output) topology file; -o (output) Run file; -p input topology file
grompp -f ${name}preEMsteep.mdp -po ${name}_EMsteep.mdp -c ${name}allHion.gro -p ${name}allHion.top -pp ${name}allHion_out.top -o ${name}allHion.tpr >& 11output.grompp_emsd2
if [ "$?" -ne "0" ] ; then
	cat 11output.grompp_emsd2
	echo grompp exited with nonzero error code
	exit 1
fi
#grompp -f minim.mdp -po ${name}_EMsteep.mdp -c ${name}allHion.gro -p ${name}allHion.top -pp ${name}allHion_out.top -o ${name}allHion.tpr >& 11output.grompp_emsd2
# steepest descent EM
# mdrun performs the molecular dynamics simulations.
# -s run input file; -o full precision trajectory output file; -c structure output file; -e energy output file; -g log file; -v "load imbalance display on screen" ?
mdrun -s ${name}allHion.tpr -o ${name}_EMsteep.trr -c ${name}_EMsteep.gro -e ${name}_EMsteep.edr -g ${name}_EMsteep.log -v >& 12output.mdrun_emsd
#if [ "$?" -ne "0" ] ; then
#	cat 12output.mdrun_emsd
#	echo mdrun exited with nonzero error code
#	exit 1
#fi
# preprocessing for EM with conjugate gradient
grompp -f ${name}preEMcg.mdp -po ${name}_EMcg.mdp -c ${name}_EMsteep.gro -p ${name}allHion.top -pp ${name}preEMsteep.top -o ${name}_EMcg.tpr >& 13output.grompp_emcg
if [ "$?" -ne "0" ] ; then
	cat 13output.grompp_emcg
	echo grompp exited with nonzero error code
	exit 1
fi
# conjugate gradient EM
mdrun -s ${name}_EMcg.tpr -o ${name}_EMcg.trr -c ${name}_EMcg.gro -e ${name}_EMcg.edr -g ${name}_EMcg.log -v >& 14output.mdrun_emcg
#if [ "$?" -ne "0" ] ; then
#	cat 14output.mdrun_emcg
#	echo mdrun exited with nonzero error code
#	exit 1
#fi


# MOLECULAR DYNAMICS
# preprecessing for PRMD
grompp -f ${name}prePRMD.mdp -po ${name}_PRMD.mdp -c ${name}_EMcg.gro -p ${name}allHion.top -pp ${name}prePRMD.top -o ${name}_PRMD.tpr -maxwarn 1 >& 15output.grommp_prmd
if [ "$?" -ne "0" ] ; then
	cat 14output.mdrun_emcg
	echo grompp exited with nonzero error code
	exit 1
fi
# position restrained and distance constrained MD
mdrun -s ${name}_PRMD.tpr -o ${name}_PRMD.trr -c ${name}_PRMD.gro -e ${name}_PRMD.edr -g ${name}_PRMD.log -v -maxh 1 >& 16output.mdrun_prmd
#if [ "$?" -ne "0" ] ; then
#	cat 16output.mdrun_prmd
#	echo mdrun exited with nonzero error code
#	exit 1
#fi
# For some reason sometimes no gro file is being generated by mdrun, so the
# following uses editconf to convert the tpr to a gro file for further use, but
# only if no gro file is detected.
if [ -f ${name}_PRMD.gro ] ; then
  echo ${name}_PRMD.gro exists already
else
  editconf -f ${name}_PRMD.tpr -o ${name}_PRMD.gro &> 16b_output.editconf
fi
# preprocessing for MD
grompp -f ${name}preMD.mdp -po ${name}_MD.mdp -c ${name}_PRMD.gro -p ${name}allHion.top -pp ${name}preMD.top -o ${name}_MD.tpr -maxwarn 1 >& 17output.grompp_md
if [ "$?" -ne "0" ] ; then
	cat 17output.grompp_md
	echo grompp exited with nonzero error code
	exit 1
fi
# MD
mdrun -s ${name}_MD.tpr -o ${name}_MD.trr -c ${name}_MD.gro -e ${name}_MD.edr -g ${name}_MD.log -v >& 18output.mdrun_md
#if [ "$?" -ne "0" ] ; then
#	cat 18output.mdrun_md
#	echo mdrun exited with nonzero error code
#	exit 1
#fi

# plot of protein's potential energy
echo 10 | g_energy -f ${name}_MD.edr -o energy.xvg >& 19output.g_energy
# generation of pdb
editconf -f ${name}_MD.gro -o ${name}_MD.pdb -bt cubic -d 0.6 >& 20output.editconf_pdb
#g_rmsf -f ${name}_MD.trr -s ${name}_MD.tpr -b 20 -e 40 -ox ${name}_MD_avg.pdb >& 21output.g_rmsf




