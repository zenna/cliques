#include <iostream>
#include <cliques/algorithms/bottleneck.h>
#include <cliques/helpers/helpers.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/structures/disjointset.h>
#include <lemon/maps.h>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#include <vector>

#include <cliques/helpers/graphhelpers.h>


//TODO
//REMOVE NEIGHBOURS TO SELF
int main() {
    lemon::SmartGraph orange_graph;

    for (int i=0;i<8;++i) {
        orange_graph.addNode();
    }

    orange_graph.addEdge(orange_graph.nodeFromId(0),orange_graph.nodeFromId(1));
    orange_graph.addEdge(orange_graph.nodeFromId(0),orange_graph.nodeFromId(2));
    orange_graph.addEdge(orange_graph.nodeFromId(1),orange_graph.nodeFromId(2));
    orange_graph.addEdge(orange_graph.nodeFromId(1),orange_graph.nodeFromId(3));
    orange_graph.addEdge(orange_graph.nodeFromId(3),orange_graph.nodeFromId(4));
    orange_graph.addEdge(orange_graph.nodeFromId(3),orange_graph.nodeFromId(5));
    orange_graph.addEdge(orange_graph.nodeFromId(4),orange_graph.nodeFromId(5));

    orange_graph.addEdge(orange_graph.nodeFromId(6),orange_graph.nodeFromId(0));
    orange_graph.addEdge(orange_graph.nodeFromId(6),orange_graph.nodeFromId(2));
    orange_graph.addEdge(orange_graph.nodeFromId(6),orange_graph.nodeFromId(1));

    orange_graph.addEdge(orange_graph.nodeFromId(7),orange_graph.nodeFromId(3));
    orange_graph.addEdge(orange_graph.nodeFromId(7),orange_graph.nodeFromId(5));
    orange_graph.addEdge(orange_graph.nodeFromId(7),orange_graph.nodeFromId(4));


    lemon::SmartGraph::Node source = orange_graph.nodeFromId(1);
	lemon::SmartGraph::Node target = orange_graph.nodeFromId(2);
	lemon::SmartGraph::EdgeMap<bool> filter(orange_graph, true);

    lemon::IterableValueMap<lemon::SmartGraph, lemon::SmartGraph::Edge, float> weights(orange_graph);

    float w = 0.5;
    for (lemon::SmartGraph::EdgeIt itr(orange_graph); itr != lemon::INVALID; ++itr) {
		weights.set(itr,w);
		w  = w + 1.0;
    }

	std::cout << "bottle neck is " << clq::find_bottleneck(orange_graph,weights,source,target,filter) << std::endl;

    return 0;
};
