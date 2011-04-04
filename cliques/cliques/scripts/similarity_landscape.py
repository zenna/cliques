"""
Generate Random Graph
Find q-measure for each graph
Find VI_matrix
Find q_dist matrix (pairwise distances in q.measure)
"""

#TODO:
# save partitions
# evaluate modularity
# Save modularity
# get positions
    
import cliques.measures.modularity as modularity
import cliques.measures.stability as stability
import networkx as nx
import cliques.ext.partition.partition_prune_fast as part
from cliques.classes.graph import Graph
from cliques.config import UNIPROJ_HOME
import pylab as P

while True:
    #G = nx.random_graphs.fast_gnp_random_graph(10, .3)
    #G = nx.read_edgelist('%s/data/ladder.pairs' % UNIPROJ_HOME)
    #G = nx.read_edgelist('%s/data/renaud_16.pairs' % UNIPROJ_HOME)
    #G = nx.complete_graph(6)
    G = nx.path_graph(7)
    
   
    #G = nx.dorogovtsev_goltsev_mendes_graph(1)
    if nx.algorithms.components.connected.is_connected(G) is True:
        break

G.__class__ = Graph

partitions = part.find_all_partitions(G.edges_as_sets())
q_list = stability.evaluate_all_partitions(G, partitions, 3.2)

import math

print "Calculating Q_matrix"
q_matrix = []
for q_1 in q_list:
    q_row = []
    for q_2 in q_list:
        q_dist = abs(q_1 - q_2) / 2
        q_dist = q_dist * math.log(len(G))
        q_row.append(q_dist)
    q_matrix.append(q_row)

print "Calculating vi-matrix"
from cliques.distances import variation_of_information
vi_matrix = variation_of_information.variation_of_information(partitions,len(G))

import numpy as np
np_q_matrix = np.array(q_matrix)
np_vi_matrix = np.array(vi_matrix)

#Prevent vi_matrix[i][j] != vi_matrix[j][i] (happens due to computuation precision)
for i in range(len(np_vi_matrix)):
    for j in range(len(np_vi_matrix)):
        np_vi_matrix[i][j] = np_vi_matrix[j][i]

sim_matrix = np_q_matrix - np_vi_matrix
sim_matrix = sim_matrix/math.log(len(G))

f_norm = np.linalg.norm(sim_matrix)
sim_matrix = abs(sim_matrix)

print "Finding Positions"
from pymatlab.matlab import MatlabSession
session = MatlabSession()
session.putvalue('vi_matrix',np_vi_matrix) 
session.run("path(path,'%s/cliques/matlab')" % UNIPROJ_HOME)
session.run("path(path,'%s/cliques/ext/modularity_visualization_v1.0.0/')" % UNIPROJ_HOME)
session.run('positions = cca(vi_matrix,2,100000)')

for i in [0.0001,0.7,0.8,1.0,10,1000000]:
    q_list = stability.evaluate_all_partitions(G, partitions, i)
    print "Calculating Q_matrix"
    q_matrix = []
    for q_1 in q_list:
        q_row = []
        for q_2 in q_list:
            q_dist = abs(q_1 - q_2) / 2
            q_dist = q_dist * math.log(len(G))
            q_row.append(q_dist)
        q_matrix.append(q_row)
        
    session.putvalue('c',np.array(q_list,dtype='Float64')) 
    session.putvalue('z',np.array(q_list,dtype='Float64')) 
    session.run('x = positions(:,1)')
    session.run('y = positions(:,2)')
    session.run('plot_surface(x,y,z,c)')
    session.run("saveas(gcf, '/home/zenna/test/path%f.fig', 'fig')"%i)
#session.putvalue('sim_matrix',sim_matrix)
#session.run('imagesc(sim_matrix)')
from cliques.draw.draw import draw_pretty
draw_pretty(G)
P.show()
print "norm is ", f_norm
