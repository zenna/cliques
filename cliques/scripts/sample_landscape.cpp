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
	typedef boost::unordered_set<VecPartition, cliques::partition_hash,
			cliques::partition_equal> VecPartitionSet;

	int num_samples = 50000;
	lemon::SmartGraph orange_graph;
	lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
	parse_arguments(ac, av, orange_graph, weights, num_samples);

	cliques::output("Sampling Partitions Uniformly");

	double markov_time = 1.0;
	double precision = 1e-9;
	int num_steps_per_sample = 3;

	cliques::find_full_normalised_stability func(orange_graph, weights,
			precision);
	std::vector<VecPartition> test = cliques::uniform_sample<VecPartition>(orange_graph,num_samples, num_steps_per_sample);

	cliques::partitions_to_file("sampling_test.txt",test);

//	cliques::sample_metropolis(orange_graph, func, markov_time,
//			num_samples, num_steps_per_sample, sample_to_count, log_uniform);
//	cliques::output("size:", sample_to_count.size());

//	for (std::unordered_map<cliques::VectorPartition, int, cliques::partition_hash>::iterator itr = sample_to_count.begin(); itr != sample_to_count.end(); ++itr) {
//		cliques::print_partition_line(itr->first);
//		double prop_energy = func(itr->first,
//					markov_time);
//		cliques::output(prop_energy,",", double(itr->second - 1) / num_samples);
//	}

	//    // Finding Maxima
	//    cliques::output("Sampling Partitions Maxima");
	//    cliques::Logging<VecPartition> log_maxima;
	//    VecPartitionSet maxima;
	//    cliques::sample_maxima(
	//    		orange_graph,
	//    		weights,
	//    		cliques::find_weighted_linearised_stability(markov_times),
	//    		cliques::linearised_stability_gain_louvain(current_markov_time),
	//    		maxima,
	//    		sampled_partitions,
	//    		log_maxima);
	//
	////   Logging<VecPartition> log_around_maxima;
	////	cliques::output("Sampling Partitions Around Maxima");
	////	int depth = 3;
	////	cliques::bfs_around_partitions(
	////			orange_graph,
	////			weights,
	////			maxima,
	////			log_around_maxima);
	////
	//
	//    VecPartitionSet all_sampled_partitions = sampled_partitions;
	//    for (auto set_itr = maxima.begin(); set_itr != maxima.end(); ++set_itr) {
	//        all_sampled_partitions.insert(*set_itr);
	//        cliques::print_partition_list(*set_itr);
	//        cliques::output("\n");
	//    }
	//
	//    arma::vec stabs_mat(stabilities.size());
	//
	//    int jj = 0;
	//    for (double t = 0.001; t< 1000; t = t * 10) {
	//    	std::vector<double> markov_times;
	//    	markov_times.push_back(t);
	//		cliques::output("Finding stabilities");
	//		cliques::find_weighted_linearised_stability compute_quality(markov_times);
	//		std::map<int, double> stabilities;
	//		for (auto set_itr = all_sampled_partitions.begin(); set_itr != all_sampled_partitions.end(); ++set_itr) {
	//			std::vector<double> stabs;
	//			cliques::Internals internals(orange_graph, weights, *set_itr);
	//			double stability = compute_quality(internals);
	//			int i = std::distance(all_sampled_partitions.begin(), set_itr);
	//			stabilities[i] = stability;
	//		}
	//
	//		for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
	//			stabs_mat(itr->first) = itr->second;
	//		}
	//    }
	//    stabs_mat.save("stabs.mat", arma::raw_ascii);
	//
	//
	//	// Distances
	//    cliques::output("Finding distances - maxima");
	//    arma::mat D_n = cliques::find_edit_dists(all_sampled_partitions);
	//
	//    cliques::output("Finding distances - rest");
	//    //arma::mat D_l = cliques::find_edit_landmark_dists(all_sampled_partitions, maxima);
	//
	//    cliques::output("finding embedding");
	//    arma::mat L = cliques::embed_mds(D_n, 3);
	//    arma::mat L_t = arma::trans(L);
	//
	//    cliques::output("saving");
	//    L_t.save("coords-sampled.mat", arma::raw_ascii);
	//
	//    auto D_y = cliques::euclid_pairwise_dists(L_t);
	//    cliques::output("residual variance", cliques::residual_variance(D_n, D_y));
	//
	//    int i=0;
	//    int num_max = maxima.size();
	//    arma::uvec maxima_mat(num_max);
	//    for (auto itr = maxima.begin(); itr != maxima.end(); ++itr) {
	//    	auto max_itr = all_sampled_partitions.find(*itr);
	//    	maxima_mat(i) = std::distance(all_sampled_partitions.begin(), max_itr);
	//    	++i;
	//    }
	//    maxima_mat.save("maxima.mat", arma::raw_ascii);
	//
	//    cliques::output("D_n",D_n.n_rows,D_n.n_cols);
	//    //cliques::output("D_l",D_l.n_rows,D_l.n_cols);
	//    cliques::output("Maxima",maxima.size());
	//    cliques::output("All",all_sampled_partitions.size());
	//
	//    cliques::output("Finding Walk");
	////    arma::umat walk = save_walk(log_all.vec, all_sampled_partitions);
	////    walk.save("walk.mat", arma::raw_ascii);
	//
	//    arma::umat walk2 = save_walk(log_uniform.vec, all_sampled_partitions);
	//    walk2.save("walk2.mat", arma::raw_ascii);
	//
	//	for (auto itr = all_sampled_partitions.begin(); itr != all_sampled_partitions.end(); ++itr) {
	//		cliques::print_partition_line(*itr);
	//	}
	//
	//    arma::colvec stabs_mat(stabilities.size());
	//    for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
	//        stabs_mat(itr->first) = itr->second;
	//    }
	//    stabs_mat.save("stabs.mat", arma::raw_ascii);


	return 0;
}
