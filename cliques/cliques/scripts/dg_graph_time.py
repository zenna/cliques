import sys
import networkx as nx
import pylab as P
from cliques.config import UNIPROJ_HOME
from pprint import pprint
import os
from cliques.distances.variation_of_information import *

#Makes a hierharchical Graph
def make_graph(partitions,markov_time,graphs,labels,all_complexities):
    i = 0
    j = 0
    partition_map = {}
    G = nx.Graph()
    
    complexity = 0.0
    for energy in sorted(partitions.iterkeys()):
        partition = partitions[energy]
        
        new_partition_map = {}
        cmplx_v = 0.0
        super_basin_size = 0
        for part in partition:
            super_basin_size = super_basin_size + len(part)
        
        for part in partition:
            G.add_node(i,{'super_basin':part, 'energy':energy})
            #labels[i] = energy
            for node in part:
                if j == 0:
                    partition_map[node] = i
                else:
                    new_partition_map[node] = i
            i = i + 1
            prob_res = float(len(part))/float(super_basin_size)
            cmplx_v = cmplx_v + prob_res * math.log(prob_res)

        print energy, cmplx_v
        complexity = complexity - cmplx_v
        #print energy, "old",len(partition_map)
        #print energy, "new",len(new_partition_map)
        if (j != 0):
            for node,part in new_partition_map.items():
                old_part = partition_map[node]
                G.add_edge(part,old_part)
            partition_map = new_partition_map.copy()
        j = j + 1
        
    graphs[markov_time] = G.copy()
    all_complexities[markov_time] = complexity
    #graphs.append(G)
    
    #P.figure(figsize=(20,20))
    #pos=nx.graphviz_layout(G,prog='dot',args='') 
    #nx.draw(G,pos,node_size=1,alpha=0.5,node_color="blue", with_labels=True,labels = labels, font_size = 7, font_color = 'red') 
    #P.axis('equal') 
    #P.savefig('circular_tree.png') 
    #P.show() 


dg_file = sys.argv[1]
graphs = {}
stabs = {}
ncs = {}
labels = {}
all_complexities = {}

all_dgs = {}
#From file create dict {time: {energy: [[]]}
lines = open(dg_file).readlines()
time_old = -999.999
energy_old = -999.99
for line in lines:
    sl = line.split(',')
    energy = float(sl[1])
    time = float(sl[0])
    maxima = int(sl[2])
    superbasin = int(sl[3])
    if time != time_old:
        all_dgs[time] = {energy:[[]]}
        time_old = time
        energy_old = energy
    elif energy != energy_old:
        all_dgs[time][energy] = [[]]
        energy_old = energy
                              
    diff = superbasin + 1 - len(all_dgs[time][energy])
    if diff > 0:
        for i in range(diff):
            #print "size, appending", superbasin, diff
            all_dgs[time][energy].append([])
    all_dgs[time][energy][superbasin].append(maxima)

for time, dg in all_dgs.items():
    print "time is :", time
    make_graph(dg,time,graphs,labels,all_complexities)
    
list_comp = []
all_times = []
all_g_sizes = []
for time in sorted(all_complexities.iterkeys()):
    all_g_sizes.append(len(graphs[time]))
    list_comp.append(all_complexities[time])
    all_times.append(time)

import numpy as np
taus = np.loadtxt('/home/zenna/repos/uniproj/data/taus/zenna_n13.taus')

from cliques.draw.draw import draw_pretty
for time in sorted(graphs.iterkeys()):
    graph = graphs[time]
    #if 0.598 <= time <= 0.65:# len(graph) > 3:
    #if 1.47 <= time <= 1.6:# len(graph) > 3:
    if len(graph) > 6:
        print time, len(graph)
        P.figure(figsize=(4,4))
        pos=nx.graphviz_layout(graph,prog='dot',args='') 
        for node, data in graph.nodes(data=True):
            pos[node] = (pos[node][0], data['energy'])
        draw_pretty(graph,pos=pos,doAxis=True)#nx.draw(graph,pos,node_size=1,alpha=0.5,node_color="blue") 
        #P.axis()
        #P.axis('equal')

P.show();
print "TIMES"
for time in sorted(all_complexities.iterkeys()):
    print time

print "COMPLEXITIES"
for time in sorted(all_complexities.iterkeys()):
    print all_complexities[time]

