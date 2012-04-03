import sys
import networkx as nx
import pylab as P
import os
from cliques.distances.variation_of_information import *
import math

folder = sys.argv[1]
#Dict Stability:partiton
partitions = {}
for root, dirs, files in os.walk(folder):
    for name in files:
        temp_partition = [[]]
        #temp_partition.append([])
        filename = os.path.join(root, name)
        print filename
        f = open(filename)
        lines = f.readlines()
        for line in lines:
            sl = line.split()
            node = int(sl[0])
            part = int(sl[1])
            
            diff = part + 1 - len(temp_partition)
            if diff > 0:
                for i in range(diff):
                    temp_partition.append([])
            temp_partition[part].append(node)
        
        partitions[float(name)] = temp_partition

G = nx.Graph()
i = 0
j = 0
skipped = 0
partition_map = {}
labels = {}

complexity = 0.0
for stability in sorted(partitions.iterkeys()):
    
    partition = partitions[stability]
    if j == 0:
        old_partition = partition
    
    new_partition_map = {}
    cmplx_v = 0.0
    super_basin_size = 0
    for part in partition:
        super_basin_size = super_basin_size + len(part)
    
    for part in partition:
        G.add_node(i,{'super_basin':part})
        labels[i] = len(part)
        for node in part:
            if j == 0:
                partition_map[node] = i
            else:
                new_partition_map[node] = i
        i = i + 1
        prob_res = float(len(part))/float(super_basin_size)
        cmplx_v = cmplx_v + prob_res * math.log(prob_res)
    complexity = complexity - cmplx_v
    print stability, "old",len(partition_map)
    print stability, "new",len(new_partition_map)
    if (j != 0):
        for node,part in new_partition_map.items():
            old_part = partition_map[node]
            G.add_edge(part,old_part)
        partition_map = new_partition_map.copy()
        old_partition = partition
        
    # Calculate All Prs
    j = j + 1    


complexity = complexity / 2
#P.figure(figsize=(20,20))
#pos=nx.graphviz_layout(G,prog='dot',args='') 
#nx.draw(G,pos,node_size=1,alpha=0.5,node_color="blue", with_labels=True,labels = labels, font_size = 7, font_color = 'red') 
#P.axis('equal') 
#P.savefig('circular_tree.png') 
##P.show() 
#
#complexity = 0.0
#for all energy levels
#    Cl = 0.0
#    for all superbasins at this energy level
#        p_res = size_of_superbasin/size_of_all_basins
#        Cl += p_res * log p_res
#    Cl = -Cl
#    complexity + Cl
#
#complexity = complexity / (energy_max - energy_min)