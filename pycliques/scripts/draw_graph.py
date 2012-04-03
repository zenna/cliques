import networkx as nx
from sys import argv
import pylab as P

if __name__ == '__main__':
    filename = argv[1]
    G = nx.read_edgelist(filename, data=(('weight',float),))
    nx.draw(G)
    P.show()
