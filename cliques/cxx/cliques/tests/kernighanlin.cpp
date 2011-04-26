#include <iostream>

#include <vector>

#include <lemon/smart_graph.h>

#include <cliques/helpers.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/kernighan_lin.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/graphhelpers.h>

int main() {
	typedef cliques::VectorPartition partition;
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<float> weights(orange_graph);

	double stability = 0;
	cliques::read_edgelist_weighted("/home/mts09/repositories/group_repository/graph-codes/cliques/data/celegansweighted.edj",
			orange_graph, weights);

	double current_markov_time = 1.0;

	std::vector<double> markov_times, adadad;
	markov_times.push_back(current_markov_time);

	std::vector<partition> optimal_partitions;
	stability = cliques::find_optimal_partition_louvain_with_gain<partition>(orange_graph,
			weights, cliques::find_weighted_linearised_stability(markov_times),
			cliques::linearised_stability_gain_louvain(current_markov_time),
			optimal_partitions);

	partition best_partition = optimal_partitions.back();
	int length = best_partition.element_count();
	for (int i = 0; i < length; i++) {
		std::cout << i << " " << best_partition.find_set(i) << "\n";
	}
	std::cout << stability << std::endl;

	cliques::VectorPartition refined_partition(lemon::countNodes(orange_graph));
	partition bad_partition(lemon::countNodes(orange_graph));
	bad_partition.initialise_as_singletons();
	cliques::find_weighted_linearised_stability compute(markov_times);
	//std::cout << "initially: " << compute(orange_graph, bad_partition, weights, adadad) << std::endl;

	double new_stability = cliques::refine_partition_kernighan_lin(
		orange_graph,
		weights,
		cliques::find_weighted_linearised_stability(markov_times),
		cliques::linearised_stability_gain_louvain(current_markov_time),
		optimal_partitions.back(),
		refined_partition);

	std::cout << new_stability << std::endl;
//	for (int i = 0; i < length; i++) {
//		std::cout << i << " " << refined_partition.find_set(i) << "\n";
//	}

	return 0;
}
