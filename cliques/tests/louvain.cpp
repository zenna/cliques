#include <iostream>
//#include <cliques/algorithms/all_partitions.h>


#include <cliques/quality_functions/stability.h>
#include <cliques/quality_functions/stability_corr.h>
#include <cliques/quality_functions/stability_info_full.h>

#include <cliques/helpers/make_graphs.h>
#include <cliques/helpers/math.h>
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

// Save partition to file
// Write multi


int main(int argc, char *argv []) {
	//std::vector<double> matrix = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
	//auto a = clq::exp(matrix,1.0,3);
	//clq::print_2d_vector(a);
	//clq::print_collection(a);


	typedef clq::VectorPartition partition;
	lemon::SmartGraph orange_graph, exp_graph;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph), exp_graph_weights(exp_graph);

	double stability = 0;
	double current_markov_time = 1;

	//clq::make_complete_graph(orange_graph,10,weights);
	clq::read_edgelist_weighted(argv[1],
			orange_graph, weights);

	std::vector<double> null_model;
	//std::vector<double> null_model = clq::create_correlation_graph_from_graph(
	//		orange_graph, weights);
	//clq::create_mutual_information_graph_from_graph(orange_graph, weights,
//			current_markov_time);

	clq::VectorPartition singletons(lemon::countNodes(orange_graph));
	singletons.initialise_as_singletons();

	//clq::find_mutual_information_stability quality;
	//clq::mutual_information_stability_gain quality_gain;
	clq::find_linearised_normalised_stability quality(current_markov_time);
	clq::linearised_normalised_stability_gain quality_gain(current_markov_time);
	clq::graph_to_exponential_graph(orange_graph, weights, exp_graph, exp_graph_weights, current_markov_time);


//	stability = quality(orange_graph, singletons, weights, singletons,null_model);
//	std::cout << "singleton stability: " << stability << std::endl;

	clq::Logging<partition> log_louvain;
	std::vector<partition> optimal_partitions;

	clq::output("Start Louvain");

	stability = clq::find_optimal_partition_louvain<partition>(
			exp_graph, exp_graph_weights, null_model, quality, quality_gain,
			singletons, optimal_partitions, 1e-9, log_louvain,true);
	clq::partitions_to_file("optimal_partitions.mat", optimal_partitions);


	partition best_partition = optimal_partitions.back();

	int length = best_partition.element_count();
	for (int i = 0; i < length; i++) {
		std::cout << i << " " << best_partition.find_set(i) << "\n";
	}
	std::cout << stability << std::endl;

	return 0;
}
