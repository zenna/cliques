#include <iostream>
//#include <cliques/algorithms/all_partitions.h>


#include <cliques/algorithms/stability.h>
//#include <cliques/algorithms/modularity.h>

#include <cliques/structures/make_graphs.h>
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

#include <cliques/graphhelpers.h>
#include <cliques/helpers.h>
//#include <cliques/helpers/math.h>


int main() {
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
	cliques::read_edgelist_weighted("/home/mts09/repositories/group_repository/graph-codes/cliques/data/demograph_stability_n640.edj", orange_graph, weights);

	cliques::VectorPartition singletons(lemon::countNodes(orange_graph));
	singletons.initialise_as_singletons();
	cliques::find_linearised_normalised_stability quality(current_markov_time);

	stability = quality(orange_graph, singletons, weights);

	cliques::Logging<partition> log_louvain;
	std::vector<partition> optimal_partitions;

	cliques::output("Start Louvain");

	stability = cliques::find_optimal_partition_louvain<partition>(
			orange_graph, weights, quality,
			cliques::linearised_normalised_stability_gain(current_markov_time),
			1e-9, singletons, optimal_partitions, log_louvain);

	partition best_partition = optimal_partitions.back();

	int length = best_partition.element_count();
	for (int i = 0; i < length; i++) {
		std::cout << i << " " << best_partition.find_set(i) << "\n";
	}
	std::cout << stability << std::endl;

	return 0;
}
