#include <iostream>
#include <cliques/algorithms/complexity.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/helpers.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/modularity.h>

#include <cliques/structures/make_graphs.h>

#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>

#include <cliques/algorithms/louvain.h>

#include <cliques/structures/vector_partition.h>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/connectivity.h>

#include <vector>

#include <cliques/graphhelpers.h>

//TODO
//REMOVE NEIGHBOURS TO SELF
int main() {
	typedef cliques::VectorPartition partition;
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<float> weights(orange_graph);
    typedef cliques::VectorPartition VecPartition;
	double stability = 0;
	cliques::make_complete_graph(orange_graph,5, weights);

	double current_markov_time = 1.0;
	cliques::VectorPartition singletons(5);
	singletons.initialise_as_singletons();

    cliques::Logging<VecPartition> log_louvain;
	std::vector<partition> optimal_partitions;
	stability = cliques::find_optimal_partition_louvain_with_gain<partition>(orange_graph,
			weights, cliques::find_weighted_linearised_stability(current_markov_time),
			cliques::linearised_stability_gain_louvain(current_markov_time),singletons,
			optimal_partitions, log_louvain);

	partition best_partition = optimal_partitions.back();
	int length = best_partition.element_count();
	for (int i = 0; i < length; i++) {
		std::cout << i << " " << best_partition.find_set(i) << "\n";
	}
	std::cout << stability << std::endl;

	return 0;
}
