#include <iostream>
//#include <cliques/algorithms/all_partitions.h>


#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/stability_corr.h>

#include <cliques/helpers/make_graphs.h>
//#include <cliques/drawing/colour_maps.h>
//
#include <cliques/algorithms/louvain.h>

#include <cliques/structures/vector_partition.h>

//#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
//#include <lemon/concepts/graph_components.h>
//#include <lemon/concepts/graph.h>
//#include <lemon/connectivity.h>

#include <vector>

#include <cliques/helpers/graphhelpers.h>
#include <cliques/helpers/helpers.h>
//#include <cliques/helpers/math.h>

// Save partition to file
// Write multi


int main(int argc, char *argv []) {
	//std::vector<double> matrix = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
	//auto a = cliques::exp(matrix,1.0,3);
	//cliques::print_2d_vector(a);
	//cliques::print_collection(a);


	typedef cliques::VectorPartition partition;
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph);

	double stability = 0;
	double current_markov_time = 1;

	//cliques::make_complete_graph(orange_graph,10,weights);
	cliques::read_edgelist_weighted(argv[1],
			orange_graph, weights);

	std::vector<double> null_model;
	//std::vector<double> null_model = cliques::create_correlation_graph_from_graph(
	//		orange_graph, weights);
	cliques::create_mutual_information_graph_from_graph(orange_graph, weights,
			current_markov_time);

	cliques::VectorPartition singletons(lemon::countNodes(orange_graph));
	singletons.initialise_as_singletons();

	cliques::find_mutual_information_stability quality;
	cliques::mutual_information_stability_gain quality_gain;

	//stability = quality(orange_graph, singletons, weights, singletons,null_model);
//	std::cout << "singleton stability: " << stability << std::endl;

	cliques::Logging<partition> log_louvain;
	std::vector<partition> optimal_partitions;

	cliques::output("Start Louvain");

	stability = cliques::find_optimal_partition_louvain<partition>(
			orange_graph, weights, null_model, quality, quality_gain,
			singletons, optimal_partitions, 1e-9, log_louvain,true);
	cliques::partitions_to_file("optimal_partitions.mat", optimal_partitions);


	partition best_partition = optimal_partitions.back();

	int length = best_partition.element_count();
	for (int i = 0; i < length; i++) {
		std::cout << i << " " << best_partition.find_set(i) << "\n";
	}
	std::cout << stability << std::endl;

	return 0;
}
