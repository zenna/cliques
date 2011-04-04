import pylab as P
import networkx as nx
#G = nx.generators.random_graphs.barabasi_albert_graph(20,2)
while True:
    temp_graph = nx.generators.random_graphs.fast_gnp_random_graph(10,.3)
    if nx.algorithms.components.connected.is_connected(temp_graph) is True:
        break
    
#G = nx.generators.random_graphs.dense_gnm_random_graph(10,15)
nx.write_edgelist(temp_graph, "barabasi.pairs", data=False)
nx.draw(temp_graph)
P.show()
