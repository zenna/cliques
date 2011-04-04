from cliques.config import UNIPROJ_HOME
import sys
#import fast_part
import networkx as nx

def find_all_connected_partitions(graph):
    graph_string = ''
    graph_string = graph_string + "%d\n" % len(graph)
    for node, degree in graph.degree().items():
        graph_string = graph_string + "%d %d\n" % (node,degree)
      
    for edge in graph.edges():
        graph_string = graph_string + "%d %d\n" % (edge[0],edge[1])
        
    return graph_string

#G = nx.read_edgelist('%s/data/renaud_16.pairs' % UNIPROJ_HOME,nodetype=int)


G = nx.classic.barbell_graph(5,4)
G = nx.complete_graph(100)
G = nx.dense_gnm_random_graph(8,9)
G = nx.moebius_kantor_graph()


G = nx.Graph()
for i in range(13):
    G.add_node(i)

G.add_edge(0,1)
G.add_edge(1,2)
G.add_edge(2,0)
G.add_edge(2,12)
G.add_edge(3,4)
G.add_edge(4,5)
G.add_edge(5,3)
G.add_edge(3,12)
G.add_edge(6,7)
G.add_edge(7,8)
G.add_edge(6,8)
G.add_edge(7,12)
G.add_edge(9,10)
G.add_edge(10,11)
G.add_edge(12,10)
G.add_edge(9,11)
G.add_edge(8,9)
#G.add_edge(7,13)
#G.add_edge(13,10)
#G.add_edge(2,14)
#G.add_edge(14,3)

#G = G = nx.dorogovtsev_goltsev_mendes_graph(4)
#G = nx.readwrite
#a = fast_part.part_wrapper()
sys.stdout.write(find_all_connected_partitions(G))
#a.find_partitions
    
#    proc = subprocess.Popen('%s/build/partition_num' % UNIPROJ_HOME, shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
#    proc.stdin.write(graph_string)
#    proc.wait()
#    output = (proc.communicate()[0])
#    return output
