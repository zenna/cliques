#include <iostream>

#include <vector>

#include <lemon/smart_graph.h>

#include <cliques/helpers.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/algorithms/kernighan_lin.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/graphhelpers.h>
#include <cliques/helpers/logger.h>

int main() {
	typedef cliques::VectorPartition partition;
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<float> weights(orange_graph);

	double stability = 0;
	cliques::read_edgelist_weighted("/home/mts09/repositories/group_repository/graph-codes/cliques/data/triangletest.edj",
			orange_graph, weights);

	double current_markov_time = 1.0;
	std::vector<double> markov_times = {1.0};

    cliques::NoLogging no_log;
	std::vector<partition> optimal_partitions;
	partition initial_partition(countNodes(orange_graph));
	initial_partition.initialise_as_singletons();
	stability = cliques::find_optimal_partition_louvain<partition>(orange_graph,
			weights, markov_times, cliques::find_linearised_normalised_stability(current_markov_time),
			cliques::linearised_normalised_stability_gain(current_markov_time),initial_partition,
			optimal_partitions, 1e-9,no_log);

	partition best_partition = optimal_partitions.back();
	std::cout << "Stab louvain " << stability << std::endl;
    cliques::print_partition_line(best_partition);
    cliques::output("End of Louvain \n");

	cliques::VectorPartition refined_partition(lemon::countNodes(orange_graph));
	partition bad_partition(lemon::countNodes(orange_graph));
	bad_partition.initialise_as_singletons();
	cliques::find_linearised_normalised_stability compute(current_markov_time);

	double new_stability = cliques::refine_partition_kernighan_lin(
		orange_graph,
		weights,
		cliques::find_linearised_normalised_stability(current_markov_time),
		cliques::linearised_normalised_stability_gain(current_markov_time),
		optimal_partitions.back(),
		refined_partition);

	std::cout << "improvement " << new_stability << std::endl;
//	for (int i = 0; i < length; i++) {
//		std::cout << i << " " << refined_partition.find_set(i) << "\n";
//	}

	return 0;
}
