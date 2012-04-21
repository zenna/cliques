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

namespace po = boost::program_options;

//template <typename V, typename S>
//arma::umat save_walk (V &vec, S &all_partitions) {
//    cliques::output("Finding Walk");
//    std::vector<int> walk_steps;
//
//    for (auto vec_itr = vec.begin(); vec_itr != vec.end(); ++vec_itr) {
//    	auto what = all_partitions.find(*vec_itr);
//    	if (what != all_partitions.end()) {
//    		int distance = std::distance(all_partitions.begin(), what);
//    		walk_steps.push_back(distance);
//    	}
//    }
//
//    arma::umat walk_mat(walk_steps.size(),2);
//    walk_mat.zeros();
//    for (unsigned int i = 0; i < walk_steps.size(); ++i) {
//    	walk_mat(i,0) = walk_steps[i];
//    }
//
//    return walk_mat;
//}

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

	cliques::output("making graph");
	if (vm.count("graph")) {
		std::string filename = vm["graph"].as<std::string> ();
		cliques::read_edgelist_weighted(filename, graph, weights);
	} else {
//		cliques::make_complete_graph(graph, 8, weights);
		//      cliques::make_ring_graph(graph, 12, weights);
		      cliques::make_complete_graph(graph, 6, weights);
	}
}

int main(int ac, char* av[]) {
	typedef cliques::VectorPartition VecPartition;
	typedef std::unordered_set<VecPartition, cliques::partition_hash,
			cliques::partition_equal> VecPartitionSet;

	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
	parse_arguments(ac, av, orange_graph, weights, num_samples);

	cliques::output("Sampling Partitions Uniformly");
	int num_samples = 50000;
	double markov_time = 1.0;
	double precision = 1e-9;
	int num_steps_per_sample = 3;

	cliques::find_full_normalised_stability func(orange_graph, weights,
			precision);
	std::vector<VecPartition> test = cliques::uniform_sample<VecPartition>(orange_graph,num_samples, num_steps_per_sample);

	cliques::output("Sampling basins");
	int num_samples_per_sample = 1000;
	vecPart initial_partition(lemon::countNodes(orange_graph));
	initial_partition.initialise_as_singletons();
	cliques::find_full_normalised_stability quality_func(orange_graph, weights, precision);

	// Use a lambda to close over orange_graph so that the function passed to stochastic_climb only takes one param
	// - the partition, but has necessary access to orange_graph

	auto unary_stochastic_monotic_climb = 
	[&orange_graph, &quality_func, &markov_time] (vecPart initial_partition) -> vecPart {
	    return cliques::stochastic_monotonic_climb
	        <vecPart, partition_set>
	        (initial_partition,
	        [&orange_graph] (vecPart partition) -> partition_set {
	            return cliques::find_neighbours(orange_graph, partition);
	        },
	        cliques::Direction::ASCENT,
	        [&quality_func, &markov_time] (vecPart partition) -> double {
	            return quality_func(partition, markov_time);
	        });
	});

	std::map<int, std::map<int, double> > basin_to_config_to_prob = 
		cliques::sample_probabalistic_basins(
			test,
			unary_stochastic_monotic_climb,
			num_samples_per_sample,
			x);


	cliques::partitions_to_file("sampling_test.txt",test);

	return 0;
}
