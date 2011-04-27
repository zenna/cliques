#include <iostream>

#include <lemon/smart_graph.h>

#include <vector>

#include <cliques/algorithms/hsg.h>
#include <cliques/graphhelpers.h>
#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/space.h>
#include <cliques/structures/vector_partition.h>

int main() {
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

	typedef cliques::VectorPartition VecPartition;
	cliques::read_edgelist_weighted(
			"/home/zenna/repos/graph-codes/cliques/data/triangletest.edj", orange_graph,
			weights);

	boost::unordered_set<VecPartition, cliques::partition_hash,
			cliques::partition_equal> all_partitions;
	cliques::find_connected_partitions(orange_graph, all_partitions);

	lemon::SmartGraph space;
	auto map = cliques::create_space(orange_graph, all_partitions, space);

	std::vector<float> stabilities;
	std::vector<double> markov_times = {1.0};
	//find_weighted_linearised_stability(markov_times);

	for (lemon::SmartGraph::EdgeIt itr(space); itr != lemon::INVALID; ++itr) {
	    stabilities.push_back(space.id(itr) % 9);
	    //need map
	    //and means to compute stability
	}

    cliques::draw_graph canvas(space);
    canvas.add_edge_map(cliques::make_energy_edge_colour_map(stabilities));
    canvas.draw("spaces.png");

	return 0;
}
