#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers.h>
#include <cliques/helpers/math.h>
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

template<typename G, typename M, typename Ps>
po::variables_map parse_arguments(int ac, char *av[], G &graph, M &weights,
		Ps &optimal_partitions, int &num_samples, int &num_dim,
		std::string &filename_prefix, bool &find_partitions,
		bool &create_space, bool &find_stabilities, bool &find_distances,
		bool &do_embedding, bool &find_basins, int &num_timesteps,
		double &start_time, double &end_time, bool merge_moveset) {
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()("help", "produce help message")("graph,G",
			po::value<std::string>(), "input graphs")("partitions_file,H",
			po::value<std::string>(), "partition file for landscape")("intermediate_graphs,I",
					po::value<std::string>(), "intermediate_graphs")(
			"num-samples,S", po::value<int>(), "number of samples")(
			"dimensions,D", po::value<int>(), "number of dimensions")(
			"prefix,x", po::value<std::string>(), "filename prefix")(
			"find_partitions,p", "Find Partitions")("create_space,s",
			"Create Space")("merge_moveset,mm", "moveset")(
			"find_stabilities,r", "Find Stabilities")("find_distances,d",
			"Find Distances")("do_embedding,e", "Do Embedding")(
			"find_basins,b", "Find Basins")("time_steps", po::value<int>(),
			"Number of time steps")("start_time", po::value<double>(),
			"start time")("end_time", po::value<double>(), "end time");

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
		cliques::output("Loading Graph");
		std::string filename = vm["graph"].as<std::string> ();
		cliques::read_edgelist_weighted(filename, graph, weights);
	} else {
		cliques::output("making default graph graph");
		//cliques::make_path_graph(graph, 7, weights);
		//      cliques::make_ring_graph(graph, 12, weights);
		cliques::make_complete_graph(graph, 8, weights);
	}

	if (vm.count("partitions_file")) {
		std::string filename = vm["partitions_file"].as<std::string> ();
		read_partitions_file(optimal_partitions, filename);
	} else {
		std::cout << "Need partitiosn file" << std::endl;
		exit(1);
	}

	return vm;
}

int main(int ac, char* av[]) {
	typedef cliques::VectorPartition VecPartition;
	typedef boost::unordered_set<VecPartition, cliques::partition_hash,
			cliques::partition_equal> VecPartitionSet;

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

	cliques::graph_to_edgelist_file(filename_prefix + "_0_graph_edgelist.edj",
			orange_graph);

	VecPartitionSet all_partitions;
	if (find_partitions) {
		cliques::output("Finding Connected Partitions");
		cliques::NoLogging no_logging;
		cliques::find_connected_partitions(orange_graph, all_partitions,
				no_logging);
		cliques::output("complete size:", all_partitions.size());

		std::ofstream vector_file;
		vector_file.open(filename_prefix + "_0_partitions.mat");
		for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
			int length = itr->element_count();
			for (int i = 0; i < length; i++) {
				vector_file << itr->find_set(i) << " ";
			}
			vector_file << std::endl;
		}
	}

	lemon::SmartGraph space;
	lemon::SmartGraph::EdgeMap<float> space_weights(space);
	boost::bimap<boost::bimaps::unordered_set_of<cliques::VectorPartition,
			cliques::partition_hash, cliques::partition_equal>,
			boost::bimaps::set_of<lemon::SmartGraph::Node> > map;
	if (create_space) {
		cliques::output("Creating space graph");
		map = cliques::create_space(orange_graph, all_partitions, space);
		lemon::SmartGraph::EdgeMap<float> space_weights(space);
		cliques::make_weights_from_edges(space, space_weights);
		cliques::output("number of nodes", lemon::countNodes(space));
		cliques::output("number of edges", lemon::countEdges(space));
	}

	std::vector<double> markov_times;
	std::vector<std::vector<double> > all_stabilities;
	if (find_stabs) {
		cliques::output("Finding stabilities");
		std::ofstream stabs_file;
		stabs_file.open(filename_prefix + "_0_energy.mat");
		cliques::output(start_time, end_time, num_timesteps);
		markov_times = cliques::create_exponential_markov_times(start_time,
				end_time, num_timesteps);
		cliques::output(markov_times.size(), "time steps");
		double precision = 1e-9;
		cliques::find_full_normalised_stability func(orange_graph, weights,
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
		cliques::graph_to_edgelist_file(
				filename_prefix + "_0_landscape_edgelist.edj", space);
	}

	arma::mat X;
	if (find_distances) {
		cliques::output("Finding distances");
		//TODO comment out appropriately..
		//auto X = cliques::find_geodesic_dists(space, landmark_nodes, space_weights);

		X = cliques::find_edit_dists(all_partitions);
		if (merge_moveset) {
			X = cliques::find_split_merge_dists(all_partitions);
			cliques::convert_dists_to_graph(space, space_weights, X, 1.0);
			cliques::graph_to_edgelist_file(
					filename_prefix + "_0_landscape_edgelist.edj", space);
		}

		X.save(filename_prefix + "_0_dists.mat", arma::raw_ascii);
	}

	arma::mat L_t;
	if (do_embedding) {
		cliques::output("finding embedding");
		auto L = cliques::embed_mds(X, num_dim);
		L_t = arma::trans(L);
		L_t.save(filename_prefix + "_0_coords.mat", arma::raw_ascii);

		auto D_y = cliques::euclid_pairwise_dists(L_t);
		cliques::output("residual variance", cliques::residual_variance(X, D_y));

	}

	// ----------------- HIERARCHY

	int size_original_graph = lemon::countNodes(orange_graph);
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

		if (cliques::read_edgelist_weighted(current_hierarchy_filename,
				tangerine_graph, weights) != true) {
			break;
		}

		cliques::graph_to_edgelist_file(filename_prefix + "_" + hierarchy_level +  + "_graph_edgelist.edj",
					tangerine_graph);

		cliques::output("Building landscape at level: ", i);

		VecPartitionSet all_partitions;
		if (find_partitions) {
			cliques::output("Finding Connected Partitions");
			cliques::NoLogging no_logging;
			cliques::find_connected_partitions(tangerine_graph, all_partitions,
					no_logging);
			cliques::output("complete size:", all_partitions.size());
			cliques::partitions_to_file(filename_prefix + "_" + hierarchy_level + "_partitions.mat",
					all_partitions);
		}

		lemon::SmartGraph space;
		lemon::SmartGraph::EdgeMap<float> space_weights(space);
		boost::bimap<boost::bimaps::unordered_set_of<cliques::VectorPartition,
				cliques::partition_hash, cliques::partition_equal>,
				boost::bimaps::set_of<lemon::SmartGraph::Node> > map2;
		if (create_space) {
			cliques::output("Creating space graph");
			map2 = cliques::create_space(tangerine_graph, all_partitions, space);
			lemon::SmartGraph::EdgeMap<float> space_weights(space);
			cliques::make_weights_from_edges(space, space_weights);
			cliques::output("number of nodes", lemon::countNodes(space));
			cliques::output("number of edges", lemon::countEdges(space));
		}

		std::vector<double> markov_times;
		std::vector<std::vector<double> > all_stabilities;
		if (find_stabs) {
			cliques::output("Finding stabilities");
			std::ofstream stabs_file;
			stabs_file.open(
					filename_prefix + "_" + hierarchy_level
							+ "_energy.mat");
			cliques::output(start_time, end_time, num_timesteps);
			markov_times = cliques::create_exponential_markov_times(start_time,
					end_time, num_timesteps);
			cliques::output(markov_times.size(), "time steps");
			double precision = 1e-9;
			cliques::find_full_normalised_stability func(tangerine_graph, weights,
					precision);
			for (unsigned int i = 0; i < markov_times.size(); ++i) {
				std::vector<double> stabilities;
				stabs_file << markov_times[i] << " ";
				for (auto itr = all_partitions.begin(); itr
						!= all_partitions.end(); ++itr) {
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
			cliques::graph_to_edgelist_file(
					filename_prefix + "_" + hierarchy_level
							+ "_landscape_edgelist.edj", space);
		}

		cliques::output("Mapping down to original level");
		arma::mat new_L_t(all_partitions.size(), num_dim);
		int h=0;
		for (auto coarse_partition = all_partitions.begin(); coarse_partition
				!= all_partitions.end(); ++coarse_partition) {
			VecPartition flat_partition(size_original_graph);
			cliques::VectorPartition &optimal_partition = optimal_partitions[i];

			// Find which partition in original landscape new partition refers to
			for (int j = 1; j < size_original_graph; j++) {
				int orig_set = optimal_partition.find_set(j);
				int coarse_set = coarse_partition->find_set(orig_set);
				flat_partition.add_node_to_set(j, coarse_set);
			}
			flat_partition.normalise_ids();
			// Find its position
			lemon::SmartGraph::Node node_in_orig_landscape = map.left.at(flat_partition);
			int id_in_orig_landscape = orange_graph.id(node_in_orig_landscape);
			new_L_t.row(h) = L_t.row(id_in_orig_landscape);
			h += 1;
		}
		new_L_t.save(filename_prefix + "_" + hierarchy_level
							+ "_coords.mat");
		i += 1;
	}

	find_basins = false;
	if (find_basins) {
		cliques::output("Finding Probabilistic Basins");
		std::vector<std::map<int, std::map<int, double>>>all_basins;
		int j = 0;
		for (auto stabilities = all_stabilities.begin(); stabilities
				!= all_stabilities.end(); ++stabilities) {
			auto basins = cliques::compute_probabalistic_basins_new(space,
					*stabilities);
			cliques::output("time", markov_times[j], "num_basins", basins.size());
			all_basins.push_back(basins);
			++j;
		}

		cliques::basins_to_file(filename_prefix + "_greedy_basins.mat", all_basins,
				markov_times);
	}

	bool find_klin_basins = false;
	if (find_klin_basins) {
		cliques::output("Finding Kerninghan Lin Basins");
		std::vector<std::map<int, std::map<int, double>>>all_basins;

		//        for (auto map_itr = map.left.begin(); map_itr != map.left.end(); ++map_itr) {
		//            cliques::print_partition_line(*map_itr);
		//        }


		for (auto time = markov_times.begin(); time != markov_times.end(); ++time) {
			lemon::SmartGraph exp_graph;
			lemon::SmartGraph::EdgeMap<double> exp_graph_weights(exp_graph);
			cliques::graph_to_exponential_graph(orange_graph, weights,exp_graph, exp_graph_weights, *time);
			//TODO map is not initialised for hasse distances.
			auto basins = cliques::compute_kernighan_lin_basins(orange_graph,weights,
					cliques::find_linearised_normalised_stability(*time),
					cliques::linearised_normalised_stability_gain(*time),
					map,
					space,
					*time,
					all_partitions);
			all_basins.push_back(basins);
		}
		cliques::basins_to_file(filename_prefix + "_klin_basins.mat", all_basins,
				markov_times);
	}
	return 0;
}

