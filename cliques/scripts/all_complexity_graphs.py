import sys
import networkx as nx
#import pylab as P
from cliques.config import UNIPROJ_HOME
import pprint
import os
from cliques.distances.variation_of_information import *
from cliques.measures.hierarchy import *
import pylab as P

#Makes a hierharchical Graph
def make_graph(root,directory,graphs,labels):
    partitions = {}
    for root, dirs, files in os.walk("%s%s" % (root,directory)):
        for name in files:
            temp_partition = [[]]
            #temp_partition.append([])
            filename = os.path.join(root, name)
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
    partition_map = {}
    for time in sorted(partitions.iterkeys()):
        partition = partitions[time]
        if j == 0:
            old_partition = partition
        
        new_partition_map = {}
        vi = vi_two_parts([old_partition,partition],100)
        if ((len(partition) < len(old_partition) and vi_two_parts([old_partition,partition],100) > math.log(100) * 0.3) ) or j == 0:
            print len(partition), len(old_partition)
            for part in partition:
                G.add_node(i)
                labels[i] = time
                for node in part:
                    if j == 0:
                        partition_map[node] = i
                    else:
                        new_partition_map[node] = i
                i = i + 1
            if (j != 0):
                for node,part in new_partition_map.items():
                    old_part = partition_map[node]
                    G.add_edge(part,old_part)
                partition_map = new_partition_map.copy()
                old_partition = partition
        j = j + 1

    p_index = int(directory[5:])
    #p_index = int(directory[6:])
    print directory, p_index, len(G)
    graphs[p_index] = G.copy()
    #graphs.append(G)
    
    #P.figure(figsize=(20,20))
    #pos=nx.graphviz_layout(G,prog='dot',args='') 
    #nx.draw(G,pos,node_size=1,alpha=0.5,node_color="blue", with_labels=True,labels = labels, font_size = 7, font_color = 'red') 
    #P.axis('equal') 
    #P.savefig('circular_tree.png') 
    #P.show() 

#Curves over all time stab
import numpy as np
def load_stab_data(root,file,stabs,ncs):
    B = np.loadtxt('%s%s' %(root,file))
    A = B[92:194]
    index = int(file[5:-7])
    #index = int(file[6:-7])
    stabs[index] = np.sum(A,0)
    ncs[index] = B[:,4]

root_folder = sys.argv[1]
G = {}
stabs = {}
ncs = {}
labels = {}
for i in range(1,33):
    filename = "COSNV%d.stdout" % i   
    load_stab_data(root_folder,filename,stabs,ncs)

for root, dirs, files in os.walk(root_folder):
    for directory in dirs:
        make_graph(root,directory,G,labels)
        

sizes = []
lengths = []
all_stabs = []

#pos=nx.graphviz_layout(G[id],prog='dot',args='')
pos= nx.circular_layout(G[1])
pos = nx.shell_layout(G[1], nlist=None, dim=2, scale=1)
#pos = spectral_layout(G[1], dim=2, weighted=True, scale=1)
x=18

import numpy

def drawHSG(x):
    y = len(G[x]) - 1
    node_groups = {}
    group_positions  = [0.0, 71.0, 144.0, 216.0, 287.0]
    pos=nx.spring_layout(G[x])
    pos=nx.graphviz_layout(G[x],prog='twopi',root=y)
    print pos[y]
    P.show()
    for node_id, pos1 in pos.items():
        #print pos[y]
        #print type(pos1)
        dist = numpy.linalg.norm(numpy.array(pos[y])-numpy.array(pos1))
        g_i = 0
        for group_position in group_positions:
            if abs(dist - group_position) < 10:
                node_groups[node_id] = g_i
            g_i = g_i + 1
    print node_groups
    edge_colours = []
    for edge in G[x].edges():
        print node_groups[edge[0]], node_groups[edge[1]]
        if node_groups[edge[0]] == node_groups[edge[1]]:
            edge_colours.append('r')
        else:
            edge_colours.append('black')
    nx.draw(G[x],pos,node_size=15,alpha=0.5,node_color="blue",edge_color=edge_colours,with_labels=False, font_color='r', font_size=10)

def drawHSGold(x):
    pos=nx.graphviz_layout(G[x],prog='dot',args='')
    nx.draw(G[x],pos,node_size=1,alpha=0.5,node_color="blue",with_labels=False, labels=labels, font_color='r', font_size=10) 
    P.axis('equal')
    P.show()

    
for id in sorted(G.iterkeys()):
    all_stabs.append(stabs[id][4])
    #sizes.append(len(G[id]))
    #lengths.append(nx.average_shortest_path_length(G[id]))
    #if id % 3 == 0:
    P.figure(figsize=(4,4))
    pos=nx.graphviz_layout(G[id],prog='dot',args='')
    for node in G[id].nodes():
        pos[node] = (pos[node][0], math.log(labels[node])*500)
    nx.draw(G[id],pos,node_size=15,alpha=0.5,node_color="blue",with_labels=False, labels=labels, font_color='r', font_size=10) 
    P.axis('equal')
    P.title(id)
    #if id > 4:
    #    break
    #P.savefig('%s/thesis/plots/hgraph%d.png' %(UNIPROJ_HOME,id)) 
    
#P.show()
#P.plot(sizes)
#P.plot(lengths)
all_nums = []
c = []
P.figure()
all_stabs = []
times = np.loadtxt('/home/zenna/repos/uniproj/data/taus/zenna_n13.taus')
for id in sorted(ncs):
    all_stabs.append(stabs[id][1]/358.0)
    all_nums.append(stabs[id][4]/358.0)
    if id <= 19:
        c.append('red')
    else:
        c.append('blue')
    P.semilogx(times,ncs[id], color = c[id-1])
    
#P.scatter(all_stabs,all_nums,color=c,label=['adada','eaea'])
#P.plot(range(1,33),all_stabs)
P.show()



noi = 31
for noi in sorted(G.iterkeys()):
    pos=nx.graphviz_layout(G[noi],prog='dot',args='')
    for node in G[noi].nodes():
        pos[node] = (math.log(labels[node]), (pos[node][0]))
    #P.figure(figsize=(7,7))
    #nx.draw(G[noi],pos,node_size=1,alpha=0.5,node_color="blue",with_labels=False)
    #P.savefig("/home/zenna/Dropbox/uniproj/brain%d.eps"%)
    


    
noi = 26
pos=nx.graphviz_layout(G[noi],prog='dot',args='')
for node in G[noi].nodes():
    pos[node] = (math.log(labels[node]), (pos[node][0]))


P.subplot(311)
P.axis('on')
P.axis((math.log(0.008), math.log(10000), -1000.0, 8000.0))
nx.draw(G[noi],pos,node_size=1,alpha=0.5,node_color="blue",with_labels=False)
P.axis('on')
P.axis((math.log(0.008), math.log(10000), -1000.0, 8000.0))

P.subplot(312)
P.axis('on')
#louvain2 = np.loadtxt('/home/zenna/repos/uniproj/data/brain_new_stabilities/COSNV%d.stdout' % noi)
louvain2 = np.loadtxt('/home/zenna/repos/uniproj/data/random_stabilities/RCOSNV%d.stdout' % noi)
P.semilogx(louvain2[:,0], louvain2[:,11])

P.subplot(313)
P.axis('on')
P.semilogx(louvain2[:,0], louvain2[:,4])
P.show();


#Robustness Calculations
variances = []
num_cs = []
weighted_v = []
shortest_paths = []
bcs = []
for id in sorted(G.iterkeys()):
    print "finding robustness", id
    #variances.append(vulnerability_variance(G[id]))
    #num_cs.append(connected_components_variance(G[id]))
    #weighted_v.append(connected_components_size_variance(G[id]))
    #shortest_paths.append(nx.average_shortest_path_length(G[id]))
    bcs.append(nx.betweenness_centrality(G[id]))
    
