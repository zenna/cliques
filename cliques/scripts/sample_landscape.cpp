#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers/make_graphs.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/helpers/logger.h>
#include <cliques/landscapes/sample_from_landscape.h>
#include <cliques/optimisation/hill_climbing.h>
#include <cliques/landscapes/non_linear_dim_reduction.h>
#include <cliques/landscapes/landscape_space.h>
#include <cliques/quality_functions/stability_full.h>
#include "cliques/landscapes/basins.h"


namespace po = boost::program_options;

template<typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights, int &num_samples) {
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()("help", "produce help message")("graph,G",
			po::value<std::string>(), "input graph")("num-samples,S",
			po::value<int>(), "number of samples");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		//return 1; - exit
	}

	if (vm.count("num-samples")) {
		num_samples = vm["num-samples"].as<int> ();
	}

	clq::output("making graph");
	if (vm.count("graph")) {
		std::string filename = vm["graph"].as<std::string> ();
		clq::read_edgelist_weighted(filename, graph, weights);
	} else {
//		clq::make_complete_graph(graph, 8, weights);
		//      clq::make_ring_graph(graph, 12, weights);
		      clq::make_complete_graph(graph, 6, weights);
	}
}

int main(int ac, char* av[]) {
	typedef typename clq::VectorPartition VecPart;
	typedef typename std::unordered_set<VecPart, clq::partition_hash, clq::partition_equal> partition_set;

	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
	int num_samples = 100;
	double start_time = 0.1, end_time = 2, num_timesteps = 2;
	parse_arguments(ac, av, orange_graph, weights, num_samples);
	std::string filename_prefix = "test";

	// Samples
	clq::output("Sampling Partitions Uniformly");
	double precision = 1e-9;
	int num_steps_per_sample = 3;
	
	auto all_partitions = clq::uniform_sample<VecPart>(orange_graph,num_samples, num_steps_per_sample);

	// Basins - must come first as may extend all_partitions
	clq::output("Sampling basins");
	clq::output(start_time, end_time, num_timesteps);
	// std::vector<double> markov_times = clq::create_exponential_markov_times(start_time,
	// 		end_time, num_timesteps);
	std::vector<double> markov_times = {1.0};
	int num_samples_per_sample = 100;
	std::mt19937 prng_engine;
	std::vector<std::map<int, std::map<int, double>>>all_basins;
	clq::output(markov_times.size(), "time steps");
	clq::find_full_normalised_stability func(orange_graph, weights,
			precision);

	// Create new bains for every time point
	for (double markov_time : markov_times) {

		// Create a new function every time to for new markov time
		// FIXME: THere is a more efficient way to do it
		// This function is basically a wrapper to the optimisation function
		// which is passed to the basin sampler
		// A wrapper is required because the basin_sampler expects an optimisation
		// function that takes one argument (the initial configuration)
		auto unary_stochastic_monotic_climb = 
		[&orange_graph, &func, &markov_time, &prng_engine]
		(VecPart initial_partition) -> VecPart {
		    return clq::stochastic_monotonic_climb
		        <VecPart, std::vector<VecPart>>
		        (initial_partition,
		        [&orange_graph] (VecPart partition) -> std::vector<VecPart> {
		            return clq::find_neighbours_vec(orange_graph, partition);
		        },
		        clq::Direction::ASCENT,
		        [&func, &markov_time] (VecPart partition) -> double {
		            return func(partition, markov_time);
		        },
		        prng_engine);
		};
		
		// This is an id function to give identity to a partition (position in container) for file export purposes
		auto get_partition_id = 
		[&all_partitions]
		(VecPart partition) -> int {
			auto find_result = std::find(all_partitions.begin(),all_partitions.end(),partition);
			if (find_result == all_partitions.end()) {
				return -1;
			}
			else {
				return std::distance(all_partitions.begin(), find_result);
			}
		};

		clq::output("Really sampling now");
		auto basins = clq::sample_probabalistic_basins
			<VecPart, std::vector<VecPart>>
			(all_partitions, unary_stochastic_monotic_climb, num_samples_per_sample, get_partition_id);

		clq::output("time", markov_time, "num_basins", basins.size());
		all_basins.push_back(basins);
	}
	clq::basins_to_file(filename_prefix + "_greedy_basins.bsn", all_basins,
			markov_times);
	
	// Must come after basin finding, since it can modify all_partitions
	clq::partitions_to_file(filename_prefix + "_partitions.mat",all_partitions);
	clq::graph_to_edgelist_file(filename_prefix + "_graph_edgelist.edj",
			orange_graph);

	// Stabilities
	std::vector<std::vector<double> > all_stabilities;
	clq::output("Finding stabilities");
	std::ofstream stabs_file;
	stabs_file.open(filename_prefix + "_energy.mat");
	for (unsigned int i = 0; i < markov_times.size(); ++i) {
		std::vector<double> stabilities;
		stabs_file << markov_times[i] << " ";
		for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
			double stability = func(*itr, markov_times[i]);
			stabilities.push_back(stability);
			stabs_file << stability << " ";
		}

		all_stabilities.push_back(stabilities);
		if (i + 1 != markov_times.size()) {
			stabs_file << std::endl;
		}
	}
	stabs_file.close();

	// Distances
	arma::mat X;
	clq::output("Finding distances");
	X = clq::find_edit_dists(all_partitions);
	X.save(filename_prefix + "_dists.mat", arma::raw_ascii);

	// Graph
	lemon::SmartGraph space;
	lemon::SmartGraph::EdgeMap<double> space_weights(space);
	clq::convert_dists_to_graph(space,space_weights,X, 1.0);
    clq::graph_to_edgelist_file(filename_prefix + "_landscape_edgelist.edj", space);
    clq::output("number of nodes", lemon::countNodes(space));
    clq::output("number of edges", lemon::countEdges(space));

	// Embedding
	clq::output("finding embedding");
	auto L = clq::embed_mds(X, 2);
	arma::mat L_t = arma::trans(L);
	L_t.save(filename_prefix + "_coords.mat", arma::raw_ascii);

	auto D_y = clq::euclid_pairwise_dists(L_t);
	clq::output("residual variance", clq::residual_variance(X, D_y));

	return 0;
}

