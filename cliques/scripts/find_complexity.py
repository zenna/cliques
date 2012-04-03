import cliques
import cliques.ext.partition.partition_prune_fast as part
from pprint import pprint
import networkx as nx
from pymatlab.matlab import MatlabSession
from time import time
all_samples = []

class MyDict(dict):
    def __missing__(self,key):
        self[key] = rv = []
        return rv

#num_samples = 100
#x = range(2,11)
#data = {}
#for num_nodes in x:
#    for num_edges in range(num_nodes-1,1+num_nodes*(num_nodes-1)/2):
#        print "starting", (num_nodes,num_edges)
#        data[(num_nodes,num_edges)] = []
#        for j in range(num_samples):
#            while True:
#                temp_graph = nx.generators.random_graphs.dense_gnm_random_graph(num_nodes,num_edges)
#                
#                if nx.algorithms.components.connected.is_connected(temp_graph) is True:
#                    break
#            
#            sample_graph = list_graph()
#            sample_graph.load_from_edges(temp_graph.edges())
#            sample_graph.find_m()
#            clusterings = part.find_all_partitions(sample_graph.graph)
#            data[(num_nodes,num_edges)].append(len(clusterings))
#
##What are the curves I want, (n,n-1) for all n
##Option 1, curve for all increments of m up to max (complete graph),
##Option 2 m = some percentage of n
#
##option 1:
#
y_data = MyDict()
for num_nodes in x:
    for num_edges in range(num_nodes-1,1+num_nodes*(num_nodes-1)/2):
        mean = np.array(data[(num_nodes,num_edges)]).mean()
        y_data[num_nodes-num_edges].append(mean)
                
session = MatlabSession()

ys = []
xs = []
for curve in y_data.values():    
    x = range(11-len(curve),11)
    x = np.array(x,dtype='Float64')
    xs.append(min(x))
    ys.append(min(curve))
    session.putvalue('x',x)
    session.putvalue('y',curve)
    session.run('semilogy(x,y)')
    session.run('hold on')  

timings = {}
clusterings_per_sec = []
for j in range (1):
    timings = {}
    clusterings_per_sec = []
    for num_nodes in range(13,15):
        #temp_graph = nx.generators.random_graphs.dense_gnm_random_graph(num_nodes,num_edges)
        while True:
            temp_graph = nx.generators.random_graphs.fast_gnp_random_graph(num_nodes,.3)
            if nx.algorithms.components.connected.is_connected(temp_graph) is True:
                break
    
        #temp_graph = nx.complete_graph(num_nodes)
        sample_graph = list_graph()
        sample_graph.load_from_edges(temp_graph.edges())
        sample_graph.find_m()
        print "starting"
        start = time()
        clusterings = part.find_all_partitions(sample_graph.graph)
        elapsed = time() - start
        timings[len(clusterings)] = elapsed
        clusterings_per_sec.append(float(len(clusterings))/float(elapsed))
        print "done %s" % num_nodes
    
print clusterings_per_sec