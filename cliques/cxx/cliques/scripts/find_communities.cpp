#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/all_communities.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/algorithms/space.h>
#include <cliques/algorithms/maxima.h>
#include <cliques/algorithms/kernighan_lin.h>
#include <cliques/algorithms/internals/internals.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/algorithms/aglob.h>

namespace po = boost::program_options;

template<typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights,
		int &num_samples, int &num_dim, std::string &filename_prefix,
		bool &find_partitions, bool &create_space, bool &find_stabilities,
		bool &find_distances, bool &do_embedding, bool &find_basins) {
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()("help", "produce help message")("graph,G",
			po::value<std::string>(), "input graph")("num-samples,S",
			po::value<int>(), "number of samples")("dimensions,D",
			po::value<int>(), "number of dimensions")("prefix,x",
			po::value<std::string>(), "filename prefix")("find_partitions,p",
			"Find Partitions")("create_space,s", "Create Space")(
			"find_stabilities,r", "Find Stabilities")("find_distances,d",
			"Find Distances")("do_embedding,e", "Do Embedding")(
			"find_basins,b", "Find Basins");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	find_partitions = vm.count("find_partitions") ? true : false;
	find_distances = vm.count("find_distances") ? true : false;

	if (vm.count("find_distances")) {
		find_distances = true;
		find_partitions = true;
	}

	if (vm.count("create_space")) {
		create_space = true;
		find_distances = true;
		find_partitions = true;
	}

	if (vm.count("find_stabilities")) {
		find_stabilities = true;
		find_partitions = true;
	}

	if (vm.count("find_basins")) {
		find_basins = true;
		create_space = true;
		find_stabilities = true;
		find_partitions = true;
		find_distances = true;
	}

	if (vm.count("do_embedding")) {
		do_embedding = true;
		find_distances = true;
		find_partitions = true;
	}

	cliques::output("bools", find_partitions, create_space, find_stabilities, find_distances, do_embedding, find_basins);

	// Do everything if nothing specified
	if (!(find_partitions || create_space || find_stabilities || find_distances
			|| do_embedding || find_basins)) {
		find_partitions = create_space = find_stabilities = find_distances
				= do_embedding = find_basins = true;
		cliques::output("F", find_partitions);
	}

	if (vm.count("help")) {
		std::cout << desc << "\n";
	}

	if (vm.count("num-samples")) {
		num_samples = vm["num-samples"].as<int> ();
	}

	if (vm.count("dimensions")) {
		num_dim = vm["dimensions"].as<int> ();
	}

	if (vm.count("prefix")) {
		filename_prefix = vm["prefix"].as<std::string> ();
	} else {
		filename_prefix = "out";
	}

	if (vm.count("graph")) {
		cliques::output("Loading Graph");
		std::string filename = vm["graph"].as<std::string> ();
		cliques::read_edgelist_weighted(filename, graph, weights);
	} else {
		cliques::output("making default graph graph");
		//cliques::make_path_graph(graph, 7, weights);
		//      cliques::make_ring_graph(graph, 12, weights);
		cliques::make_complete_graph(graph, 8, weights);
	}
}

int main(int ac, char* av[]) {
	typedef cliques::VectorPartition VecPartition;
	typedef boost::unordered_set<VecPartition, cliques::partition_hash,
			cliques::partition_equal> VecPartitionSet;

	lemon::SmartGraph orange_graph;
	std::string filename_prefix;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
	int num_samples = 100000;
	int num_dim = 3;
	bool find_partitions, create_space, find_stabs, find_basins,
			find_distances, do_embedding;
	find_partitions = create_space = find_stabs = find_basins =
	            find_distances = do_embedding = false;
	parse_arguments(ac, av, orange_graph, weights, num_samples, num_dim,
			filename_prefix, find_partitions, create_space, find_stabs,
			find_distances, do_embedding, find_basins);

	cliques::graph_to_edgelist_file(filename_prefix + "_graph_edgelist.edj",
			orange_graph);

//	std::vector<int> comm1 = {0,1};
//	std::vector<int> comm2 = {0,1};
//	std::vector<int> comm3 = {2};
//	std::vector<int> buffer(4,0);
//
//	auto d1 = cliques::find_community_dist(orange_graph, weights, comm1, comm2, buffer);
//	auto d2 = cliques::find_community_dist(orange_graph, weights, comm1, comm3, buffer);
//
//	cliques::output("dists", d1, d2);



//	markov_times = cliques::create_exponential_markov_times(0.00001, 500, 500);
	cliques::find_full_normalised_stability func(orange_graph, weights,1e-9);
//	double stability = func(p, 1, markov_times[i]);
//	std::vector<std::vector<double> > all_stabilities;
//	auto all_maxima = cliques::find_optimal_communities_huxley(orange_graph, weights, func, 1.0, all_stabilities);
//	cliques::output("maximum");
//	for (auto maximum = all_maxima.begin(); maximum != all_maxima.end(); ++maximum) {
//		cliques::print_collection(*maximum);
//	}

    return 0;
}
