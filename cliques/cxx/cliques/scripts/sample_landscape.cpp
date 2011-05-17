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

int main() {

    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    typedef cliques::VectorPartition VecPartition;

//    cliques::output("making graph");
//    cliques::make_path_graph(orange_graph, 16);
    cliques::make_complete_graph(orange_graph, 8);
//    cliques::read_edgelist_weighted(
//            "/home/zenna/repos/graph-codes/cliques/data/graphs/renaud_n12.edj",
//            orange_graph, weights);

    typedef cliques::VectorPartition VecPartition;

    cliques::output("Sampling Partitions");
    boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> sampled_partitions;
    cliques::sample_space_uniform(orange_graph,1000, 10, sampled_partitions);


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
    auto L = cliques::find_embedding_mds_smacof(X, 3);
    arma::mat L_t = arma::trans(L);

    cliques::output("saving");
    L_t.save("coords.mat", arma::raw_ascii);

    auto D_y = cliques::euclid_pairwise_dists(L_t);
    cliques::output("residual variance", cliques::residual_variance(X, D_y));
    return 0;
}
