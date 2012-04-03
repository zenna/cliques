import community
import networkx as nx
import pylab as plt
from numpy import *

#better with karate_graph() as defined in networkx example.
#erdos renyi don't have true community structure
#G = nx.erdos_renyi_graph(30, 0.05)
G = nx.read_edgelist('/home/mts09/workspace/MeanField/data/networks/test.txt')
#first compute the best partition
partition = community.best_partition(G)

#drawing
size = float(len(set(partition.values())))

# pos = nx.spring_layout(G)
print size
pos = {}
for node in partition.keys():
    ind = str(node)	
    xcoord = 10*cos(partition[node]*2*pi/size) + 1*cos(float(node)*2*pi/size) #random.random()
    ycoord = 10*sin(partition[node]*2*pi/size) + 1*sin(float(node)*2*pi/size)#random.random()
    pos[ind] = array([xcoord,ycoord])

count = 0.
for com in set(partition.values()) :
    count = count + 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 20,
                                node_color = str(count / size))


nx.draw_networkx_edges(G,pos, alpha=0.5)
plt.show()

