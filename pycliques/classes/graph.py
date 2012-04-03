import networkx as nx

class Graph(nx.Graph):
    def __init__(self,filename=None):
        self.m = 0
        self.graph = {}
        if filename is not None:
            self.load(filename)

    def find_m(self):
        for (node, connected_nodes) in self.graph.items():
            self.m = self.m + len(connected_nodes)

    #Adjacency matrix - Assumes Undirected
    def A(self,one,two):
        if two in self.edges_as_sets[one]:
            return 1
        else:
            return 0

    #Degree of a node
    def k(self,node):
        return len(self.edges_as_sets[node])
    
    def edges_as_sets(self):
        self.edges_as_sets = {}
        edges = self.edges()
        for edge in edges:
            u = int(edge[0]) #networkx importer casts as string, but partitioner needs ints
            v = int(edge[1])
            if u not in self.edges_as_sets:
                self.edges_as_sets[u] = set()
            if v not in self.edges_as_sets:
                self.edges_as_sets[v] = set()
            self.edges_as_sets[u].add(v)
            self.edges_as_sets[v].add(u)
        
        return self.edges_as_sets
        
    def load(self,filename):
        try:
            for line in file(filename):
                line = line.split()
                u = int(line[0])
                v = int(line[1])
                if u not in self.graph:
                    self.graph[u] = set()
                if v not in self.graph:
                    self.graph[v] = set()
                self.graph[u].add(v)
                self.graph[v].add(u)
        except IOError:
            print "Could not find filename %s:" % filename