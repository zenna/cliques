#include <iostream>
#include <algorithms/complexity.h>
#include <algorithms/partitions.h>
#include <helpers.h>
#include <algorithms/stability.h>
#include <algorithms/modularity.h>

#include <drawing/draw.h>
#include <drawing/colour_maps.h>

#include <algorithms/module.h>
//TODO get rid of this
#include <structures/disjointset.h>

#include <cliques/structures/vector_partition.h>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/connectivity.h>

#include <vector>

#include <graphhelpers.h>

//TODO
//REMOVE NEIGHBOURS TO SELF
int main() {
	typedef cliques::VectorPartition partition;
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<float> weights(orange_graph);

	cliques::read_edgelist_weighted("/home/mts09/celegansweighted.edj",
			orange_graph, weights);

	double current_markov_time = 0.8;

	std::vector<partition> optimal_partitions;
	cliques::find_optimal_partition_louvain_with_gain<partition>(orange_graph,
			weights,
			cliques::linearised_stability_louvain(current_markov_time),
			cliques::linearised_stability_gain_louvain(current_markov_time),
			optimal_partitions);

	partition best_partition = optimal_partitions.back();
//	int length = best_partition.element_count();
//	for (int i = 0; i < length; i++) {
//		std::cout << i << " " << best_partition.find_set(i) << "\n";
//	}

	/*	cliques::print_partition(optimal_partitions.back());

	 //Drawing
	 float start = -10;
	 std::vector<float> energies;
	 for (lemon::SmartGraph::EdgeIt e(orange_graph); e != lemon::INVALID; ++e) {
	 energies.push_back(start);
	 start += 12.0;
	 std::cout << orange_graph.id(e) << std::endl;
	 }

	 cliques::draw_graph canvas(orange_graph);
	 canvas.add_node_map(cliques::make_partition_colour_map<
	 cliques::DisjointSetForest<int> >(best_partition));
	 canvas.add_edge_map(cliques::make_energy_edge_colour_map(energies));
	 canvas.draw("test_louvain_out");
	 */
	return 0;
}
