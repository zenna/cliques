#!/bin/bash

# script for doing all the fancy landscape operations in one go

if [ -z $1 ] ; then
echo "you need to provide one argument to identify paths"
exit 0
fi

if [ $1 = "M" ] ; then
# Michael's DIR
DIR_REP="$HOME/repositories/cliques"
DIR_BIN="$HOME/PhD/ext_programmes_and_lib/cliques"
fi

if [ $1 = "Z" ] ; then
# Zenna's DIR
DIR_REP="$HOME/repos/cliques"
DIR_BIN="$HOME/builds/cliques/debug"
fi

# check if argument 2 is empty if yes use default test graph otherwise use graph provided
if [ -z $2 ] ; then
GRAPH="barbell_n8.edj"
else 
GRAPH=$2
fi

echo "Run: $DIR_BIN/tests/test_louvain $DIR_REP/data/graphs/$GRAPH"
echo `$DIR_BIN/tests/test_louvain $DIR_REP/data/graphs/$GRAPH`
echo "Louvain algorithm finished"

echo `$DIR_BIN/scripts/create_multilandscape -G $DIR_REP/data/graphs/$GRAPH -I ./intermediate_graphs -H ./optimal_partitions.mat \n`
GRAPH=${GRAPH%.edj}
GRAPH=${GRAPH##*/}
echo "create multilandscape finished"


echo `python $DIR_REP/pycliques/scripts/create_json.py -m -x ./out -o $DIR_REP/jscliques/data/$GRAPH.json \n`
echo "OUTPUT FILE WRITTEN TO $DIR_REP/jscliques/data/$GRAPH.json"

if [ -n $3  ] ; then
	if [ "$3" = "R" ] ; then
		echo `mkdir output_$GRAPH`
		echo `mv out_* output_$GRAPH/` 
		echo `mv intermediate_graphs* output_$GRAPH/` 
		echo `mv optimal_partitions.mat output_$GRAPH/`
	fi
fi 
