import networkx as nx
import numpy as np

def efficiency(graph):
    return 1.0/nx.average_shortest_path_length(graph)

def vulnerability_variance(graph):
    efficiencies = []
    for node in graph.nodes():
        temp_graph = graph.copy()
        eff_old = efficiency(temp_graph)
        temp_graph.remove_node(node)
        eff_new = efficiency(temp_graph)
        relative_eff = (eff_old - eff_new) / eff_old
        efficiencies.append(relative_eff)
        
    return np.var(efficiencies)

def connected_components_variance(graph):
    num_connected_components = []
    for node in graph.nodes():
        temp_graph = graph.copy()
        temp_graph.remove_node(node)
        num_connected_components.append(nx.number_connected_components(temp_graph))
        
    return np.max(num_connected_components)

def connected_components_size_variance(graph):
    weighted_vulnerability = 0.0
    for node in graph.nodes():
        temp_graph = graph.copy()
        temp_graph.remove_node(node)
        comps = nx.connected_components(temp_graph)
        if len(comps) > 1:
            comp_lengths = []
            for comp in comps:
                comp_lengths.append(len(comp))
                
            print "num,var:", len(comps), np.var(comp_lengths)
            #print weighted_vulnerability
            weighted_vulnerability = weighted_vulnerability + np.exp(-np.var(comp_lengths)/10000.0)
        
    return weighted_vulnerability

#
#G = nx.Graph()
#B = nx.Graph()
#
#for i in range(8):
#    G.add_node(i)
#    B.add_node(i)
#    
#G.add_edge(0,1)
#G.add_edge(1,3)
#G.add_edge(1,4)
#G.add_edge(1,5)
#G.add_edge(0,2)
#G.add_edge(2,6)
#G.add_edge(2,7)
#G.add_edge(2,8)
#
#B.add_edge(0,1)
#B.add_edge(1,3)
#B.add_edge(1,4)
#B.add_edge(1,5)
#B.add_edge(0,2)
#B.add_edge(2,6)
#B.add_edge(2,7)
#B.add_edge(2,8)
#B.add_edge(1,2)
##
##D = nx.empty_graph(8)
##print vulnerability_variance(G)
##print vulnerability_variance(B)
#print connected_components_size_variance(G)
#print connected_components_size_variance(B)