#include <iostream>
#include <algorithms/partitions.h>
#include <helpers.h>
#include <algorithms/stability.h>

#include <drawing/draw.h>
#include <drawing/colour_maps.h>

#include <lemon/list_graph.h>

#include <vector>

#include <graphhelpers.h>


//TODO
//REMOVE NEIGHBOURS TO SELF
int main() {

	// Create Graph
    lemon::ListGraph orange_graph;

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

	// Find all partitions
    cliques::umap partition_map;
    find_connected_partitions(orange_graph, partition_map);

    //Create landscape

    //Find stabilities


	// Create edge weights from graph -> edge_weights (bottleneck)
	std::vector<float> current_markov_time;
    current_markov_time.push_back(0.2);
    current_markov_time.push_back(1.0);
    current_markov_time.push_back(1000);
    std::vector<float> stabilities;

    cliques::draw_graph canvas(orange_graph);
    //canvas.draw(cliques::make_energy_edge_colour_map(edge_weights));

    return 0;
};
