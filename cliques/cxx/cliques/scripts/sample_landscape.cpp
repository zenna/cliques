#include <vector>
#include <iostream>
#include <algorithm>

#include <lemon/smart_graph.h>

#include <cliques/graphhelpers.h>
#include <cliques/algorithms/stability.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>
#include <cliques/algorithms/space.h>
#include <cliques/algorithms/maxima.h>
#include <cliques/structures/vector_partition.h>


struct NoLogging {
    template <typename P>
    void log(const P &) {}
};

template <typename P>
struct Logging {
    std::vector<P> vec;

    //template <typename P>
    void log(const P &p) {
        vec.push_back(p);
    }
};

int main() {
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    typedef cliques::VectorPartition VecPartition;

    cliques::output("making graph");
    cliques::make_path_graph(orange_graph, 4, weights);
//    cliques::make_complete_graph(orange_graph, 4, weights);
//    cliques::read_edgelist_weighted(
//            "/home/zenna/repos/graph-codes/cliques/data/graphs/renaud_n12.edj",
//            orange_graph, weights);

    Logging<VecPartition> log_uniform;
    cliques::output("Sampling Partitions Uniformly");
    boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> sampled_partitions;
    cliques::sample_uniform_metropolis(orange_graph,10000, 10, sampled_partitions, log_uniform);

//    Logging<VecPartition> log_louvain;
//    cliques::output("Sampling Partitions Louvain");
//	std::vector<partition> optimal_partitions;
	double current_markov_time = 1.0;
	std::vector<double> markov_times = {1.0};
//    cliques::find_optimal_partition_louvain_with_gain<partition>(orange_graph,
//    			weights, cliques::find_weighted_linearised_stability(markov_times),
//    			cliques::linearised_stability_gain_louvain(current_markov_time),
//    			optimal_partitions, log_lovauin);
//
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

    Logging<VecPartition> log_maxima;
    boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> maxima;
    cliques::output("Sampling Partitions Maxima");
    cliques::sample_maxima(
    		orange_graph,
    		weights,
    		cliques::find_weighted_linearised_stability(markov_times),
    		cliques::linearised_stability_gain_louvain(current_markov_time),
    		maxima,
    		sampled_partitions,
    		log_maxima);

//    Logging<VecPartition> log_around_maxima;
//	cliques::output("Sampling Partitions Around Maxima");
//	int depth = 3;
//	cliques::bfs_around_partitions(
//			orange_graph,
//			weights,
//			maxima,
//			log_around_maxima);

    boost::unordered_set<VecPartition, cliques::partition_hash,
                cliques::partition_equal> all_sampled_partitions = sampled_partitions;

    for (auto set_itr = maxima.begin(); set_itr != maxima.end(); ++set_itr) {
        all_sampled_partitions.insert(*set_itr);
        cliques::print_partition_list(*set_itr);
        cliques::output("\n");
    }

    cliques::output("Finding stabilities");
    cliques::find_weighted_linearised_stability compute_quality(markov_times);
    std::map<int, double> stabilities;
    for (auto set_itr = all_sampled_partitions.begin(); set_itr != all_sampled_partitions.end(); ++set_itr) {
        std::vector<double> stabs;
        cliques::Internals internals(orange_graph, weights, *set_itr);
        cliques::print_partition_list(*set_itr);
        cliques::output("comm_w_in");
        for (int i =0;i < internals.comm_w_in.size(); ++i) {
            cliques::output(internals.comm_w_in[i]);
        }
        cliques::output("comm_w_tot");
        for (int i =0;i < internals.comm_w_tot.size(); ++i) {
            cliques::output(internals.comm_w_tot[i]);
        }
        cliques::output("node_to_w");
        for (int i =0;i < internals.node_to_w.size(); ++i) {
            cliques::output(internals.node_to_w[i]);
        }
        cliques::output("node weight to com");
        for (int i =0;i < internals.node_weight_to_communities.size(); ++i) {
            cliques::output(internals.node_weight_to_communities[i]);
        }
        cliques::output("two_m", internals.two_m);
        cliques::output("num_nodes",internals.num_nodes);
        double stability = compute_quality(internals);
        int i = std::distance(all_sampled_partitions.begin(), set_itr);
        //cliques::output(std::distance(maxima.begin(), set_itr));
        stabilities[i] = stability;
    }

    arma::colvec stabs_mat(stabilities.size());
    for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
        stabs_mat(itr->first) = itr->second;
    }
    stabs_mat.save("stabs.mat", arma::raw_ascii);

    cliques::output("Finding distances - maxima", maxima.size());
    arma::mat D_n = cliques::find_edit_dists(all_sampled_partitions);
    cliques::output("Finding distances - rest");
    arma::mat D_l = cliques::find_edit_landmark_dists(all_sampled_partitions, maxima);
    cliques::output("finding embedding");
    arma::mat L = cliques::embed_mds(D_n, 2);
    arma::mat L_t = arma::trans(L);
    L_t.save("trilaterated.mat", arma::raw_ascii);

    cliques::output("D_n",D_n.n_rows,D_n.n_cols);
    cliques::output("D_l",D_l.n_rows,D_l.n_cols);
    cliques::output("Maxima",maxima.size());
    cliques::output("All",all_sampled_partitions.size());

//    cliques::output("Finding Walk");
//    int num_steps = log.vec.size();
//    arma::umat walk(num_steps,2);
//    walk.zeros();
//    for (int i=0;i<num_steps;++i) {
//        auto what = sampled_partitions.find(log.vec[i]);
//        walk(i,0) = std::distance(all_sampled_partitions.begin(), what);
//    }
//    walk.save("walk.mat", arma::raw_ascii);



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
