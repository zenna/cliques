import sys
sys.path.append('ext/partition/')
import partition_prune_fast as part
from pprint import pprint

#Input to this file a list of cluterings
# each element will be a list, the elements of which will be lists representing clusters
# each element of which will be a node numbera 
# [ [1,2], [4,5,6] ]

# goal is to find modularity in each clustering
# need to know the degree of each node
    # can do this by calculating from 

# Currently runs clusterings code, finds all all clusterings as global variable and uses that
# Need to reorganise into modules such that can switch in/out 
    # 1. Clustering scheme, perhaps I want random clusters, or not all clusters
    # 2. What I do with the clustering(s), check modularity, check stability,
        #. Need to interface with matlab to do stability clustering

#Solutio:
#Have a module for clusterings where each class is a different clustering method
# and has a method find_clusterings which returns the clusterings

#Have another module for the analysis of clusterings
#Compare all clusterings
        
     
 

    
    
degree = part.g
clusterings = part.b
mods = []

m = 0
for (node, connected_nodes) in degree.items():
    m = m + len(connected_nodes)

#Adjacency matrix - Assumes Undirected
def A(one,two):
    if two in degree[one]:
        return 1
    else:
        return 0
    
def k(node):
    return float(len(degree[node]))
    
def find_modularity():
    global mods
    global clusterings
    global m	# total number of links -i.e. #lines in incidence list
    for clustering in clusterings:
        Q = float(0)
        for cluster in clustering:
            for i in cluster:
                for j in cluster:
                    Q = Q + float(A(i,j)) - (k(i)*k(j))/(float(m))
    
        Q = Q/(m)
        mods.append((Q, clustering))
    
    mods = sorted(mods)
    pprint(mods)
            
find_modularity()
print mods