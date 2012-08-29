#include <iostream>
#include <vector>

#include <lemon/smart_graph.h>

#include "armadillo"

#include <cliques/quality_functions/stability.h>
#include <cliques/quality_functions/stability_corr.h>
#include <cliques/quality_functions/stability_info_full.h>
#include <cliques/helpers/make_graphs.h>
#include <cliques/helpers/math.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/helpers/graphhelpers.h>
#include <cliques/helpers/helpers.h>
#include "cliques/algorithms/varofinf.h"

int main(int argc, char *argv []) {
	typedef clq::VectorPartition partition;
	lemon::SmartGraph orange_graph, exp_graph;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph), exp_graph_weights(exp_graph);

	double stability = 0;
	double current_markov_time = 5.0;

	clq::read_edgelist_weighted(argv[1],
			orange_graph, weights);

	std::vector<double> null_model;
	clq::VectorPartition singletons(lemon::countNodes(orange_graph));
	singletons.initialise_as_singletons();

	clq::find_linearised_normalised_stability quality(current_markov_time);
	clq::linearised_normalised_stability_gain quality_gain(current_markov_time);
	clq::graph_to_exponential_graph(orange_graph, weights, exp_graph, exp_graph_weights, current_markov_time);

	clq::Logging<partition> log_louvain;
	std::vector<std::vector<partition>> all_optimal_partitions;

	int num_iterations = 100;
	clq::output("Start Louvain");
	int max_levels_seen = -10, min_levels_seen = 100000;
	// Find optimal_partitions_for number of iterations
	for (int i=0;i<num_iterations;++i) {
		std::vector<partition> optimal_partitions;
		stability = clq::find_optimal_partition_louvain<partition>(
				exp_graph, exp_graph_weights, null_model, quality, quality_gain,
				singletons, optimal_partitions, 1e-9, log_louvain,true);
		all_optimal_partitions.push_back(optimal_partitions);
		clq::output("size", optimal_partitions.size());
		int size = optimal_partitions.size();
		if (size > max_levels_seen) {
			max_levels_seen = optimal_partitions.size();
		}
		if (size < min_levels_seen) {
			min_levels_seen = optimal_partitions.size();
		}
	}
	clq::partitions_to_file("optimal_partitions.mat", all_optimal_partitions[0]);

	clq::output("max_levels_seen: ",max_levels_seen, "min_levels_seen:", min_levels_seen);
	
	for (int level = 0; level < max_levels_seen; ++level) {
		arma::mat X(num_iterations,num_iterations);
		for (int i=0;i<all_optimal_partitions.size();++i) {
			for (int j=0; j<all_optimal_partitions.size();++j) {
				if (all_optimal_partitions[i].size() > level &&
					all_optimal_partitions[j].size() > level) {
					double var_of_inf = clq::find_variation_of_information(all_optimal_partitions[i][level],all_optimal_partitions[j][level]);
					X(i,j) = var_of_inf;
				}
				else {
					X(i,j) = 1.5;
				}
			}
		}
		std::stringstream level_ss;
		level_ss << level;
		X.save("vi_lovain_dists_"+ level_ss.str() +".csv", arma::raw_ascii);
	}
	return 0;
}
