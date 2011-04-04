from cliques.partition.all_connected import *
import networkx as nx

G = nx.read_edgelist('%s/data/renaud_16.pairs' % UNIPROJ_HOME,nodetype=int)
op = find_all_connected_partitions(G)
print op