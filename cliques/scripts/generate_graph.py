from cliques.config import UNIPROJ_HOME
import networkx as nx
import subprocess
from time import time
import pickle
import pylab as P

all_results = []
results = []
times = []
num_samples = 1

absolute_start = time()

for j in range(num_samples):
    results = []
    for i in range (2,18):
        graph = ''
        while True:
            #G = nx.complete_graph(i)
            p = 3.0/float(i)
            if p >= 1:
                p = 0.999
            G = nx.generators.random_graphs.fast_gnp_random_graph(i,p)
        #G = nx.read_edgelist('%s/data/renaud_16.pairs' % UNIPROJ_HOME,nodetype=int)
        #G = nx.complete_graph(3)
      
            if nx.algorithms.components.connected.is_connected(G) is True:
                break

        graph = graph + "%d\n" % len(G)
        for node, degree in G.degree().items():
            #print node, degree
            graph = graph + "%d %d\n" % (node,degree)
          
        for edge in G.edges():
            graph = graph + "%d %d\n" % (edge[0],edge[1])
            #print edge[0], edge[1]
        
        start = time()
        proc = subprocess.Popen('%s/build/partition_num' % UNIPROJ_HOME, shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.stdin.write(graph)
        proc.wait()
        elapsed = time() - start
        print j, i
        num_partitions = int(proc.communicate()[0])
        results.append(num_partitions)
        times.append([num_partitions, elapsed])
        #results.append([i, num_partitions, elapsed])
        print "%s" % (num_partitions)
        #rint "Rate is %d clusters per second" % (float(num_partitions)/float(elapsed))
    all_results.append(results)

import numpy as np
A = np.array(all_results)
#np.savetxt('num_parts.dat',A)
P.semilogy(range(2,18),np.mean(A,0))

T = np.array(times)
#np.savetxt('num_parts_times.dat',T)
P.figure()
P.loglog(T[:,0],T[:,1], 'o')

f = open('%s/data/renaud_partitions.pickle' % UNIPROJ_HOME,'w')
#pickle.dump(all_results,f)
f.close()

import pylab as P
import random
print "doing edge comparison"
num_nodes = 10
edge_results = []
all_edge_results = []

pos=nx.graphviz_layout(nx.generators.random_graphs.dense_gnm_random_graph(num_nodes,i),prog="neato")
for j in range(1):
    edge_results = []
    for i in range (num_nodes - 1, ((num_nodes * (num_nodes - 1))/2) + 1):
        graph = ''
        while True:
            G = nx.generators.random_graphs.dense_gnm_random_graph(num_nodes,i)
            
            if nx.algorithms.components.connected.is_connected(G) is True:
                break


        if i%4 == 0:
            P.figure()
            #pos=nx.graphviz_layout(G,prog="neato")
            nx.draw(G,pos,node_color='#A0CBE2',edge_color='#A0CBE2',width=4,with_labels=False)
             
        graph = graph + "%d\n" % len(G)
        for node, degree in G.degree().items():
            graph = graph + "%d %d\n" % (node,degree)
          
        for edge in G.edges():
            graph = graph + "%d %d\n" % (edge[0],edge[1])
        
        start = time()
        proc = subprocess.Popen('%s/build/partition_num' % UNIPROJ_HOME, shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        proc.stdin.write(graph)
        proc.wait()
        elapsed = time() - start
        print j, i
        num_partitions = int(proc.communicate()[0])
        edge_results.append(num_partitions)
        #results.append([i, num_partitions, elapsed])
        print "%s" % (num_partitions)
        #rint "Rate is %d clusters per second" % (float(num_partitions)/float(elapsed))
    all_edge_results.append(edge_results)

D = np.array(all_edge_results)
#np.savetxt('num_parts_edges.dat',A)
P.figure()
P.semilogy(range (num_nodes - 1, ((num_nodes * (num_nodes - 1))/2) + 1), np.mean(D,0))
P.show()
total_elapsed = time() - absolute_start
print "total time is ", total_elapsed 
#for result in results:
#  print float(result[1])/result[2]
