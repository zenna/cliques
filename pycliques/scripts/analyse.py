import sys
sys.path.append('ext/partition/')
import partition_prune_fast as part
from pprint import pprint
from measures.modularity import *
#from measures.stability import *
from graph import *
import os
UNIPROJ_HOME = os.getenv('UNIPROJ_HOME')

#Compare Partitionings
#Sort results
#interface with matlab
#

#TODO
#WANT CURVE FOR CLUSTERINGS I GIVE IT
#i.e. stability for each time step for each clustering
#Can create stability plot for different time steps
#Recreate technique used in X, compared with all clustering technique
#Show how their technique compares with a random sampling
#Create stability curve for X, does it also have the degeneracies


graph = list_graph('%s/data/test_graph_5_4.dat'%UNIPROJ_HOME)
graph.find_m()
clusterings = part.find_all_partitions(graph.graph)
mod = Modularity()
mods = mod.evaluate_all_clusterings(graph,clusterings)

#stab = Stability()
#stabs = stab.evaluate_all_clusterings(graph,clusterings)
