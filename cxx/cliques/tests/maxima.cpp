#include <iostream>

#include <algorithms/maxima.h>
#include <helpers.h>

#include <lemon/maps.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>

#include <vector>
#include <set>


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
    float *stabilities = new float[8];
    stabilities[0] = 3.0;
    stabilities[1] = 10.0;
    stabilities[2] = 9.0;
    stabilities[6] = 5.0;
    stabilities[3] = 8.0;
    stabilities[4] = 5.0;
    stabilities[5] = 8.5;
    stabilities[7] = 11.0;

    cliques::print_collection(cliques::find_maxima(orange_graph, stabilities));

    //std::vector<int> alpha;
    //cliques::print_vector(alpha);

	//std::cout << "bottle neck is " << cliques::find_bottleneck(orange_graph,weights,source,target,filter) << std::endl;

    return 0;
};
