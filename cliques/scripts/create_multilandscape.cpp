#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers/helpers.h>
#include <cliques/helpers/math.h>
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

template<typename G, typename M, typename Ps>
po::variables_map parse_arguments(int ac, char *av[], G &graph, M &weights,
		Ps &optimal_partitions, int &num_samples, int &num_dim,
		std::string &filename_prefix, bool &find_partitions,
		bool &create_space, bool &find_stabilities, bool &find_distances,
		bool &do_embedding, bool &find_basins, int &num_timesteps,
		double &start_time, double &end_time, bool merge_moveset) {
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()("help", "produce help message")
		("graph,G", po::value<std::string>(), "input graphs")
		("partitions_file,H", po::value<std::string>(), "partition file for landscape")
		("intermediate_graphs,I", po::value<std::string>(), "intermediate_graphs")
		("num-samples,S", po::value<int>(), "number of samples")
		("dimensions,D", po::value<int>(), "number of dimensions")
		("prefix,x", po::value<std::string>(), "filename prefix")
		("find_partitions,p", "Find Partitions")
		("create_space,s","Create Space")("merge_moveset,mm", "moveset")
		("find_stabilities,r", "Find Stabilities")
		("find_distances,d","Find Distances")("do_embedding,e", "Do Embedding")
		("find_basins,b", "Find Basins")
		("num_timesteps", po::value<int>(),"Number of time steps")
		("start_time", po::value<double>(),"start time")
		("end_time", po::value<double>(), "end time");

	po::variables_map vm;
	po::store(po::parse_command_line(ac, av, desc), vm);
	po::notify(vm);

	find_partitions = vm.count("find_partitions") ? true : false;
	find_distances = vm.count("find_distances") ? true : false;

	if (vm.count("time_steps")) {
		num_timesteps = vm["time_steps"].as<int> ();
	}

	if (vm.count("start_time")) {
		start_time = vm["start_time"].as<double> ();
	}

	if (vm.count("end_time")) {
		end_time = vm["end_time"].as<double> ();
	}

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

	// Do everything if nothing specified
	if (!(find_partitions || create_space || find_stabilities || find_distances
			|| do_embedding || find_basins)) {
		find_partitions = create_space = find_stabilities = find_distances
				= do_embedding = find_basins = true;
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

	if (vm.count("merge_moveset")) {
		merge_moveset = true;
		create_space = false;
	}

	if (vm.count("graph")) {
		clq::output("Loading Graph");
		std::string filename = vm["graph"].as<std::string> ();
		clq::read_edgelist_weighted_maintain_ids(filename, graph, weights);
	} else {
		clq::output("making default graph graph");
		//clq::make_path_graph(graph, 7, weights);
		//      clq::make_ring_graph(graph, 12, weights);
		clq::make_complete_graph(graph, 8, weights);
	}

	if (vm.count("partitions_file")) {
		std::string filename = vm["partitions_file"].as<std::string> ();
		read_partitions_file(optimal_partitions, filename);
	} else {
		std::cout << "Need partitions file" << std::endl;
		exit(1);
	}

	return vm;
}


template <typename G, typename W>
void do_work(G &graph,
	W &weights,
	std::string prefix,
	int num_dim,
	double start_time,
	double end_time,
	int num_timesteps,
	bool create_space,
	bool find_partitions,
	bool find_stabs,
	bool find_distances,
	bool do_embedding,
	bool find_basins) {

	typedef clq::VectorPartition VecPartition;
	typedef std::unordered_set<VecPartition, clq::partition_hash,
			clq::partition_equal> VecPartitionSet;

	clq::graph_to_edgelist_file(prefix + "_graph_edgelist.edj",
			graph);
	
	VecPartitionSet all_partitions;
	if (find_partitions) {
		clq::output("Finding Connected Partitions");
		clq::NoLogging no_logging;
		clq::find_connected_partitions(graph, all_partitions,
				no_logging);
		clq::output("complete size:", all_partitions.size());
		partitions_to_file(prefix + "_partitions.mat", all_partitions);
	}

	lemon::SmartGraph space;
	lemon::SmartGraph::EdgeMap<float> space_weights(space);
	boost::bimap<boost::bimaps::unordered_set_of<clq::VectorPartition,
			clq::partition_hash, clq::partition_equal>,
			boost::bimaps::set_of<lemon::SmartGraph::Node> > map;
	if (create_space) {
		clq::output("Creating space graph");
		map = clq::create_space(graph, all_partitions, space);
		lemon::SmartGraph::EdgeMap<float> space_weights(space);
		clq::make_weights_from_edges(space, space_weights);
		clq::output("number of nodes", lemon::countNodes(space));
		clq::output("number of edges", lemon::countEdges(space));
	}

	std::vector<double> markov_times;
	std::vector<std::vector<double> > all_stabilities;
	if (find_stabs) {
		clq::output("Finding stabilities");
		std::ofstream stabs_file;
		stabs_file.open(prefix + "_energy.mat");
		clq::output(start_time, end_time, num_timesteps);
		markov_times = clq::create_exponential_markov_times(start_time,
				end_time, num_timesteps);
		clq::output(markov_times.size(), "time steps");
		double precision = 1e-9;
		clq::find_full_normalised_stability func(graph, weights,
				precision);
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
		clq::graph_to_edgelist_file(prefix + "_landscape_edgelist.edj", space);
	}

	arma::mat X;
	if (find_distances) {
		clq::output("Finding distances");
		//TODO comment out appropriately..
		//auto X = clq::find_geodesic_dists(space, landmark_nodes, space_weights);

		X = clq::find_edit_dists(all_partitions);
		// if (merge_moveset) {
		// 	X = clq::find_split_merge_dists(all_partitions);
		// 	clq::convert_dists_to_graph(space, space_weights, X, 1.0);
		// 	clq::graph_to_edgelist_file(
		// 			prefix + "_landscape_edgelist.edj", space);
		// }

		X.save(prefix + "_dists.mat", arma::raw_ascii);
	}

	arma::mat L_t;
	if (do_embedding) {
		clq::output("finding embedding");
		auto L = clq::embed_mds(X, num_dim);
		L_t = arma::trans(L);
		L_t.save(prefix + "_coords.mat", arma::raw_ascii);
		// auto D_y = clq::euclid_pairwise_dists(L_t);
		// clq::output("residual variance", clq::residual_variance(X, D_y));
	}

	if (find_basins) {
		clq::output("Finding Probabilistic Basins");
		std::vector<std::map<int, std::map<int, double>>>all_basins;
		int j = 0;
		for (auto stabilities : all_stabilities) {
			auto basins = clq::compute_mp_probabalistic_basins(space,
					stabilities);
			clq::output("time", markov_times[j], "num_basins", basins.size());
			all_basins.push_back(basins);
			++j;
		}

		clq::basins_to_file(prefix + "_greedy_basins.bsn", all_basins,
				markov_times);
	}
}

int main(int ac, char* av[]) {
	typedef clq::VectorPartition VecPartition;
	typedef std::unordered_set<VecPartition, clq::partition_hash,
			clq::partition_equal> VecPartitionSet;

	lemon::SmartGraph orange_graph;
	std::string filename_prefix;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
	std::vector<VecPartition> optimal_partitions;
	int num_samples = 100000, num_dim = 3, num_timesteps = 500;
	double start_time = 1e-5, end_time = 500;
	bool find_partitions, create_space, find_stabs, find_basins,
			find_distances, do_embedding, merge_moveset;
	find_partitions = create_space = find_stabs = find_basins = find_distances
			= do_embedding = merge_moveset = false;

	po::variables_map vm = parse_arguments(ac, av, orange_graph, weights, optimal_partitions,
			num_samples, num_dim, filename_prefix, find_partitions,
			create_space, find_stabs, find_distances, do_embedding,
			find_basins, num_timesteps, start_time, end_time, merge_moveset);

	clq::graph_to_edgelist_file(filename_prefix + "_0_graph_edgelist.edj",
			orange_graph);

	do_work(orange_graph,
		weights,
		filename_prefix + "_0",
		num_dim,
		start_time,
		end_time,
		num_timesteps,
		create_space,
		find_partitions,
		find_stabs,
		find_distances,
		do_embedding,
		find_basins);

	std::string hierarchy_prefix = vm["intermediate_graphs"].as<std::string> ();
	int i = 1;
	while (true) {
		std::stringstream hierarchy_level_sstream;
		hierarchy_level_sstream << i;
		std::string hierarchy_level(hierarchy_level_sstream.str());

		std::string current_hierarchy_filename = hierarchy_prefix
				+ "_" + hierarchy_level;
		lemon::SmartGraph tangerine_graph;
		lemon::SmartGraph::EdgeMap<double> weights(tangerine_graph);

		if (clq::read_edgelist_weighted_maintain_ids(current_hierarchy_filename,
				tangerine_graph, weights) != true) {
			break;
		}

		do_work(tangerine_graph,
			weights,
			filename_prefix + "_" + hierarchy_level,
			num_dim,
			start_time,
			end_time,
			num_timesteps,
			create_space,
			find_partitions,
			find_stabs,
			find_distances,
			do_embedding,
			find_basins);


		i += 1;
	}	


	return 0;
}

