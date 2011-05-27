#include <vector>
#include <iostream>
#include <algorithm>

#include <lemon/smart_graph.h>

#include <cliques/algorithms/all_partitions.h>
#include <cliques/graphhelpers.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/louvain.h>

#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>
#include <cliques/algorithms/space.h>
#include <cliques/algorithms/maxima.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/helpers/logger.h>

template <typename V, typename S>
arma::umat save_walk (V &vec, S &all_partitions) {
    cliques::output("Finding Walk");
    std::vector<int> walk_steps;

    for (auto vec_itr = vec.begin(); vec_itr != vec.end(); ++vec_itr) {
    	auto what = all_partitions.find(*vec_itr);
    	if (what != all_partitions.end()) {
    		int distance = std::distance(all_partitions.begin(), what);
    		walk_steps.push_back(distance);
    	}
    }

    arma::umat walk_mat(walk_steps.size(),2);
    walk_mat.zeros();
    for (unsigned int i = 0; i < walk_steps.size(); ++i) {
    	walk_mat(i,0) = walk_steps[i];
    }

    return walk_mat;
}

int main() {
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    typedef cliques::VectorPartition VecPartition;
    typedef boost::unordered_set<VecPartition, cliques::partition_hash,
                cliques::partition_equal> VecPartitionSet;

    cliques::output("making graph");
//    cliques::make_path_graph(orange_graph, 8, weights);
//    cliques::make_ring_graph(orange_graph, 12, weights);

//    cliques::make_complete_graph(orange_graph, 6, weights);
    cliques::read_edgelist_weighted(
            "/home/zenna/repos/graph-codes/cliques/data/graphs/barbell_n8.edj",
            orange_graph, weights);

    cliques::output("Sampling Partitions Uniformly");
    cliques::Logging<VecPartition> log_uniform;
    VecPartitionSet sampled_partitions;
    cliques::sample_uniform_metropolis(orange_graph,5000, 10, sampled_partitions, log_uniform);

    cliques::output("size:",sampled_partitions.size());

//    // Sample Partitions
//    cliques::output("Finding Louvain Partitions");
//    cliques::Logging<VecPartition> log_louvain;
//    std::vector<VecPartition> optimal_partitions;
//	double current_markov_time = 0.001;
//	std::vector<double> markov_times = {0.001};
//    cliques::find_optimal_partition_louvain_with_gain<VecPartition>(orange_graph,
//    			weights, cliques::find_weighted_linearised_stability(markov_times),
//    			cliques::linearised_stability_gain_louvain(current_markov_time),
//    			optimal_partitions, log_louvain);

//    Logging<VecPartition> log_klin;
//    cliques::output("Sampling Partitions Kernighan Lin");
//    VecPartition input_partition, output_partition;
//    input_partition.initialise_as_singletons();
//    cliques::refine_partition_kernighan_lin(
//    		orange_graph,
//    		weights,
//    		cliques::find_weighted_linearised_stability(markov_times),
//    		cliques::linearised_stability_gain_louvain(current_markov_time),
//    		input_partition,
//    		output_partition,
//    		log_klin);

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
//
//	VecPartitionSet all_sampled_partitions = sampled_partitions;
////    for (auto set_itr = maxima.begin(); set_itr != maxima.end(); ++set_itr) {
////        all_sampled_partitions.insert(*set_itr);
////        cliques::print_partition_list(*set_itr);
////        cliques::output("\n");
////    }
//
////
////    arma::vec stabs_mat(stabilities.size());
////
////    int jj = 0;
////    for (double t = 0.001; t< 1000; t = t * 10) {
////    	std::vector<double> markov_times;
////    	markov_times.push_back(t);
////		cliques::output("Finding stabilities");
////		cliques::find_weighted_linearised_stability compute_quality(markov_times);
////		std::map<int, double> stabilities;
////		for (auto set_itr = all_sampled_partitions.begin(); set_itr != all_sampled_partitions.end(); ++set_itr) {
////			std::vector<double> stabs;
////			cliques::Internals internals(orange_graph, weights, *set_itr);
////			double stability = compute_quality(internals);
////			int i = std::distance(all_sampled_partitions.begin(), set_itr);
////			stabilities[i] = stability;
////		}
////
////		for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
////			stabs_mat(itr->first) = itr->second;
////		}
////    }
////    stabs_mat.save("stabs.mat", arma::raw_ascii);
//
//	// Distances
//    cliques::output("Finding distances - maxima", maxima.size());
//    arma::mat D_n = cliques::find_edit_dists(all_sampled_partitions);
//    cliques::output("Finding distances - rest");
//    //arma::mat D_l = cliques::find_edit_landmark_dists(all_sampled_partitions, maxima);
//    cliques::output("finding embedding");
//    arma::mat L = cliques::embed_mds(D_n, 3);
//    arma::mat L_t = arma::trans(L);
//    L_t.save("trilaterated.mat", arma::raw_ascii);
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
////
////    arma::colvec stabs_mat(stabilities.size());
////    for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
////        stabs_mat(itr->first) = itr->second;
////    }
////    stabs_mat.save("stabs.mat", arma::raw_ascii);
//
//    cliques::output("Finding distances");
//    auto X = cliques::find_edit_dists(sampled_partitions);
//
//    cliques::output("finding embedding");
//    auto L = cliques::find_embedding_mds_smacof(X, 2);
//    arma::mat L_t = arma::trans(L);
//
//    cliques::output("saving");
//    L_t.save("coords.mat", arma::raw_ascii);
//
//    auto D_y = cliques::euclid_pairwise_dists(L_t);
//    cliques::output("residual variance", cliques::residual_variance(X, D_y));
    return 0;
}
