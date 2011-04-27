#include <iostream>

#include <lemon/smart_graph.h>

#include <vector>

#include <cliques/algorithms/hsg.h>
#include <cliques/graphhelpers.h>
#include <cliques/drawing/draw.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/space.h>
#include <cliques/structures/vector_partition.h>

int main() {
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

	typedef cliques::VectorPartition VecPartition;
	cliques::read_edgelist_weighted(
			"/home/zenna/repos/uniproj/data/triangle.edj", orange_graph,
			weights);

	boost::unordered_set<VecPartition, cliques::partition_hash,
			cliques::partition_equal> all_partitions;
	cliques::find_connected_partitions(orange_graph, all_partitions);

	lemon::SmartGraph space;
	cliques::create_space(orange_graph, all_partitions, space);

    cliques::draw_graph canvas(space);
    //canvas.add_node_map(cliques::make_partition_colour_map<cliques::DisjointSetForest<int> >(partition));
    canvas.draw("spaces.png");

	return 0;
}
