#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers/helpers.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/all_communities.h>
#include <cliques/quality_functions/stability.h>
#include <cliques/quality_functions/stability_full.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/landscapes/landscape_space.h>
#include <cliques/algorithms/kernighan_lin.h>
#include <cliques/quality_functions/internals/internals.h>
#include <cliques/landscapes/non_linear_dim_reduction.h>
#include <cliques/helpers/make_graphs.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/algorithms/aglob.h>
#include "cliques/landscapes/basins.h"

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

	clq::output("bools", find_partitions, create_space, find_stabilities, find_distances, do_embedding, find_basins);

	// Do everything if nothing specified
	if (!(find_partitions || create_space || find_stabilities || find_distances
			|| do_embedding || find_basins)) {
		find_partitions = create_space = find_stabilities = find_distances
				= do_embedding = find_basins = true;
		clq::output("F", find_partitions);
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
		clq::output("Loading Graph");
		std::string filename = vm["graph"].as<std::string> ();
		clq::read_edgelist_weighted(filename, graph, weights);
	} else {
		clq::output("making default graph graph");
		//clq::make_path_graph(graph, 7, weights);
		//      clq::make_ring_graph(graph, 12, weights);
		clq::make_complete_graph(graph, 8, weights);
	}
}

int main(int ac, char* av[]) {
	double precision = 1e-9;
	typedef clq::VectorPartition VecPartition;
	typedef std::unordered_set<VecPartition, clq::partition_hash,
			clq::partition_equal> VecPartitionSet;

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

	clq::graph_to_edgelist_file(filename_prefix + "_graph_edgelist.edj",
			orange_graph);

	std::vector<std::vector<int> > communities;
	if (find_partitions) {
		clq::output("Find Connected Communities");
		communities = clq::find_connected_communities(orange_graph);
		clq::output(communities.size(), "communities");
	    std::ofstream vector_file;
	    vector_file.open(filename_prefix + "_partitions.mat");
	    for (auto itr = communities.begin(); itr != communities.end(); ++itr) {
	        clq::VectorPartition p = clq::community_to_partition(
	                orange_graph, *itr, 0);
	        int length = p.element_count();
	        for (int i = 0; i < length; i++) {
	            vector_file << p.find_set(i) << " ";
	        }
	        vector_file << std::endl;
	    }
	}

	std::vector<double> markov_times;
	std::vector<std::vector<double> > all_stabilities;
	if (find_stabs) {
		clq::output("Finding stabilities");
		std::ofstream stabs_file;
		stabs_file.open(filename_prefix + "_energy.mat");
		clq::find_full_normalised_stability func(orange_graph, weights, precision);
		markov_times = clq::create_exponential_markov_times(0.00001, 500, 500);
		clq::output("time_steps", markov_times.size());
		for (unsigned int i = 0; i < markov_times.size(); ++i) {
			clq::output("time",markov_times[i]);
			std::vector<double> stabilities;
			stabs_file << markov_times[i] << " ";
			for (auto itr = communities.begin(); itr != communities.end(); ++itr) {
				clq::VectorPartition p = clq::community_to_partition(
						orange_graph, *itr, 0);
				double stability = func(p, 1, markov_times[i]);
				stabilities.push_back(stability);
				stabs_file << stability << " ";
			}

			all_stabilities.push_back(stabilities);
			if (i + 1 != markov_times.size()) {
				stabs_file << std::endl;
			}
		}
		stabs_file.close();
	}

	arma::mat X;
	if (find_distances) {
		clq::output("Finding distances");
		X = clq::find_community_edit_dists(orange_graph, communities);
		X.save(filename_prefix + "_dists.mat", arma::raw_ascii);
	}
//    clq::output("Finding distances Hasse");
//    clq::create_hasse_community_space(space,communities,space);
//    std::vector<lemon::SmartGraph::Node> all_nodes;
//    for (unsigned int i = 0; i < communities.size(); ++i) {
//        all_nodes.push_back(space.nodeFromId(i));
//    }
//    clq::make_weights_from_edges(space,space_weights);
//    auto X = clq::convert_graph_to_geodesic_dists(space,all_nodes,space_weights);

	lemon::SmartGraph space;
	lemon::SmartGraph::EdgeMap<double> space_weights(space);
	if (create_space) {
		clq::convert_dists_to_graph(space,space_weights,X, 1.0);
	    clq::graph_to_edgelist_file(filename_prefix + "_landscape_edgelist.edj", space);
	    clq::output("number of nodes", lemon::countNodes(space));
	    clq::output("number of edges", lemon::countEdges(space));
	}

	if (do_embedding) {
		clq::output("Finding embedding");
		auto L = clq::embed_mds(X, num_dim);
		arma::mat L_t = arma::trans(L);
		L_t.save(filename_prefix + "_coords.mat", arma::raw_ascii);

		auto D_y = clq::euclid_pairwise_dists(L_t);
		clq::output("residual variance", clq::residual_variance(X, D_y));
	}

	if (find_basins) {
		clq::output("Finding Probabalistic Community Basins");
		std::vector<std::map<int, std::map<int, double> >> all_basins;
		int j =0;
		for (auto stabilities = all_stabilities.begin(); stabilities != all_stabilities.end(); ++ stabilities) {
			auto basins = clq::compute_mp_probabalistic_basins(space, *stabilities);
			clq::output("time", markov_times[j], "num_basins",basins.size());
			all_basins.push_back(basins);
			++j;
		}
		clq::basins_to_file(filename_prefix + "_greedy_basins.bsn", all_basins, markov_times);
	}
    return 0;
}
