#include <vector>
#include <iostream>
#include <algorithm>

#include <lemon/smart_graph.h>

#include <cliques/graphhelpers.h>
#include <cliques/algorithms/stability.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>
#include <cliques/algorithms/space.h>
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
//    cliques::make_path_graph(orange_graph, 16);
    cliques::make_complete_graph(orange_graph, 7);
//    cliques::read_edgelist_weighted(
//            "/home/zenna/repos/graph-codes/cliques/data/graphs/renaud_n12.edj",
//            orange_graph, weights);

    Logging<VecPartition> log_uniform;
    cliques::output("Sampling Partitions Uniformly");
    boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> sampled_partitions;
    cliques::sample_uniform_metropolis(orange_graph,100000, 1, sampled_partitions, log_uniform);

    Logging<VecPartition> log_louvain;
    cliques::output("Sampling Partitions Louvain");
	std::vector<partition> optimal_partitions;
	double current_markov_time = 1.0;
	std::vector<double> markov_times = {1.0};
    cliques::find_optimal_partition_louvain_with_gain<partition>(orange_graph,
    			weights, cliques::find_weighted_linearised_stability(markov_times),
    			cliques::linearised_stability_gain_louvain(current_markov_time),
    			optimal_partitions, log_lovauin);

    Logging<VecPartition> log_klin;
    cliques::output("Sampling Partitions Kernighan Lin");
    VecPartition input_partition, output_partition;
    input_partition.initialise_as_singletons();
    cliques::refine_partition_kernighan_lin(
    		orange_graph,
    		weights,
    		cliques::find_weighted_linearised_stability(markov_times),
    		cliques::linearised_stability_gain_louvain(current_markov_time),
    		input_partition,
    		output_partition,
    		log_klin);

    Logging<VecPartition> log_maxima;
    cliques::output("Sampling Partitions Maxima");
    cliques::find_maxima_sampled(
    		orange_graph,
    		weights,
    		cliques::find_weighted_linearised_stability(markov_times),
    		cliques::linearised_stability_gain_louvain(current_markov_time),
    		maxima,
    		log_maxima);

    Logging<VecPartition> log_around_maxima;
	cliques::output("Sampling Partitions Around Maxima");
	int depth = 3;
	cliques::bfs_around_partitions(
			orange_graph,
			weights,
			maxima,
			log_around_maxima);

	//Set union
	// Produce landmrks from set of maxima + random sampling
	// Embed
	// Compute stabilities / save
	// Compute walks / save

    cliques::output("Finding Walk");
    int num_steps = log.vec.size();
    arma::umat walk(num_steps,2);
    walk.zeros();
    for (int i=0;i<num_steps;++i) {
        auto what = sampled_partitions.find(log.vec[i]);
        walk(i,0) = std::distance(sampled_partitions.begin(), what);
    }
    walk.save("walk.mat", arma::raw_ascii);

//    cliques::output("Finding stabilities");
//    std::vector<double> markov_times = { 1.0 };
//    cliques::find_weighted_linearised_stability func(markov_times);
//    std::map<int, double> stabilities;
//    for (lemon::SmartGraph::NodeIt itr(space); itr != lemon::INVALID; ++itr) {
//        std::vector<double> stabs;
//        VecPartition p = map.right.at(itr);
//        cliques::Internals internals(orange_graph, weights, p);
//        double stability = func(internals);
//        stabilities[orange_graph.id(itr)] = stability;
//    }
//
//    arma::colvec stabs_mat(stabilities.size());
//    for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
//        stabs_mat(itr->first) = itr->second;
//    }
//    stabs_mat.save("stabs.mat", arma::raw_ascii);

    cliques::output("Finding distances");
    auto X = cliques::find_edit_dists(sampled_partitions);

    cliques::output("finding embedding");
    auto L = cliques::find_embedding_mds_smacof(X, 2);
    arma::mat L_t = arma::trans(L);

    cliques::output("saving");
    L_t.save("coords.mat", arma::raw_ascii);

    auto D_y = cliques::euclid_pairwise_dists(L_t);
    cliques::output("residual variance", cliques::residual_variance(X, D_y));
    return 0;
}
