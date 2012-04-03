#!/bin/bash

# script for doing all the fancy landscape operations in one go

if [ -z $1 ] ; then
echo "you need to provide one argument to identify paths"
exit 0
fi

if [ $1 = "M" ] ; then
# Michael's DIR
DIR_REP="$HOME/repositories/group_repository/graph-codes/cliques"
DIR_BIN="$HOME/PhD/ext_programmes_and_lib/clique"
fi

if [ $1 = "Z" ] ; then
# Zenna's DIR
DIR_REP="$HOME/repos/graph-codes/cliques"
DIR_BIN="$HOME/cliques2"
fi

# check if argument 2 is empty if yes use default test graph otherwise use graph provided
if [ -z $2 ] ; then
GRAPH="barbell_n8.edj"
else 
GRAPH=$2
fi

echo `$DIR_BIN/tests/test_louvain $DIR_REP/data/graphs/$GRAPH`
echo `$DIR_BIN/scripts/create_multilandscape -G $DIR_REP/data/graphs/$GRAPH -I ./intermediate_graphs -H ./optimal_partitions.mat`
echo `python $DIR_REP/cliques/scripts/create_json.py -x ./out -o $DIR_REP/cliques.js/interactive/js/data/$GRAPH.json`

if [ -n $3  -a $3 = "R" ] ; then
	echo `rm out_* intermediate_graphs* optimal_partitions.mat`
fi 
