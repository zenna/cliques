#!/bin/bash

# script for doing all the fancy landscape operations in one go

# Michael's DIR
DIR_REP="$HOME/repositories/cliques"
DIR_BIN="$HOME/PhD/ext_programmes_and_lib/cliques"
DIR_WKSP="$HOME/workspace/landscapes"

# check if argument 2 is empty if yes use default test graph otherwise use graph provided
if [ -z $1 ] ; then
echo "you need to supply a graph file"
else 
GRAPH=$1

fi

echo "Run: $DIR_BIN/tests/test_louvain $DIR_WKSP/$GRAPH"
echo `$DIR_BIN/tests/test_louvain $DIR_WKSP/$GRAPH`
echo "Louvain algorithm finished"

echo `$DIR_BIN/scripts/create_multilandscape -G $DIR_WKSP/$GRAPH -I ./intermediate_graphs -H ./optimal_partitions.mat \n`
GRAPH=${GRAPH%.edj}
GRAPH=${GRAPH##*/}
echo "create multilandscape finished"


echo `mkdir output_$GRAPH`
echo `mv out_* output_$GRAPH/` 
echo `mv intermediate_graphs* output_$GRAPH/` 
echo `mv optimal_partitions.mat output_$GRAPH/`
