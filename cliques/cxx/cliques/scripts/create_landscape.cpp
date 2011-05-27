#include <vector>
#include <iostream>
#include <algorithm>

#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>

#include <cliques/graphhelpers.h>
#include <cliques/drawing/draw.h>
#include <cliques/drawing/colour_maps.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/stability.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>

#include <cliques/algorithms/space.h>
#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/structures/vector_partition.h>

// TODO: I'm not convinced maxima is working correctly
// TODO: Restrain find_neighbours to connected (or not)
// TODO: Hijack Louvain Properly

int main() {

    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    typedef cliques::VectorPartition VecPartition;
    typedef boost::unordered_set<VecPartition, cliques::partition_hash,
                cliques::partition_equal> VecPartitionSet;

    cliques::output("making graph");
    //cliques::make_path_graph(orange_graph, 3, weights);
    //cliques::make_ring_graph(orange_graph, 12, weights);
    //cliques::make_complete_graph(orange_graph, 3);
    cliques::read_edgelist_weighted(
            "/home/zenna/repos/graph-codes/cliques/data/graphs/barbell_n8.edj",
            orange_graph, weights);

    cliques::output("Finding Connected Partitions");
    cliques::Logging<VecPartition> log_all;
    VecPartitionSet all_partitions;
    int num_partitions = cliques::find_connected_partitions(orange_graph, all_partitions, log_all);

    cliques::output("complete size:",num_partitions);

//    cliques::output("Finding stabilities");
//    std::vector<double> markov_times = { 1.0 };
//    cliques::find_weighted_linearised_stability func(markov_times);
//    double stability = 0.0;
//    for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
//        VecPartition p = *itr;
//        cliques::Internals internals(orange_graph, weights, p);
//        stability = func(internals);
//        cliques::output("stability", stability);
//    }

//    std::cout << "Computing Space" << std::endl;
//    lemon::SmartGraph space;
//    auto map = cliques::create_space(orange_graph, all_partitions, space);
//    lemon::SmartGraph::EdgeMap<float> space_weights(space);
//    for (lemon::SmartGraph::EdgeIt e(space); e!= lemon::INVALID; ++e) {
//        space_weights[e] = 1.0;
//    }
//
//    cliques::output("Finding stabilities");
//    std::vector<double> markov_times = { 1.0 };
//    cliques::find_weighted_linearised_stability func(markov_times);
//    std::map<int, double> stabilities;
//    for (lemon::SmartGraph::NodeIt itr(space); itr != lemon::INVALID; ++itr) {
//        std::vector<double> stabs;
//        VecPartition p = map.right.at(itr);
//        cliques::LinearisedInternals internals(orange_graph, weights, p);
//        double stability = func(internals);
//        stabilities[orange_graph.id(itr)] = stability;
//    }
//
//    arma::colvec stabs_mat(stabilities.size());
//    for (auto itr = stabilities.begin(); itr != stabilities.end(); ++itr) {
//        stabs_mat(itr->first) = itr->second;
//    }
//    stabs_mat.save("stabs.mat", arma::raw_ascii);
//
//    cliques::output("num partitions:", num_partitions);
//    cliques::output("choosing nodes");
////    auto landmark_nodes = cliques::randomly_choose_nodes(space, 500);
//    std::vector<lemon::SmartGraph::Node> landmark_nodes;
//    for (lemon::SmartGraph::NodeIt n(space); n!= lemon::INVALID; ++n) {
//        lemon::SmartGraph::Node node = n;
//        landmark_nodes.push_back(node);
//    }
//
//    cliques::output("Finding distances");
//    //auto X = cliques::find_geodesic_dists(space, landmark_nodes, space_weights);
//    //auto Y = cliques::find_edit_dists(space, landmark_nodes, map);
//    auto X = cliques::find_edit_dists(all_partitions);
//
//    cliques::output("finding embedding");
//    auto L = cliques::embed_mds(X, 2);
//    arma::mat L_t = arma::trans(L);
//
//    //auto D_y = cliques::euclid_pairwise_dists(L_t);
//    //cliques::output("num partitions:", num_partitions);
//    //cliques::output("residual variance", cliques::residual_variance(X, D_y));
//
//    L_t.save("coords.mat", arma::raw_ascii);
    /*std::vector<VecPartition> optimal_partitions;
    stability
            = cliques::find_optimal_partition_louvain_with_gain<VecPartition>(
                    space, space_weights,
                    cliques::find_weighted_linearised_stability(markov_times),
                    cliques::linearised_stability_gain_louvain(
                            current_markov_time), optimal_partitions);

    VecPartition best_partition = optimal_partitions.back();
    best_partition.normalise_ids();
    std::cout << "num_sets: " << *std::max_element(
            best_partition.partition_vector.begin(),
            best_partition.partition_vector.end());

    std::vector<float> stabilities;
     std::vector<double> markov_times = { 1.0 };
     cliques::find_weighted_linearised_stability func(markov_times);

     std::cout << "computing stability" << std::endl;
     for (lemon::SmartGraph::NodeIt itr(space); itr != lemon::INVALID; ++itr) {
     std::vector<double> stabs;
     VecPartition p = map.right.at(itr);
     cliques::Internals internals(orange_graph, weights, p);
     double stability = func(internals);
     stabilities.push_back(stability);
     std::cout << stability << std::endl;
     }

     std::cout << "num_nodes" << lemon::countNodes(space) << std::endl;
     std::cout << "num_edges" << lemon::countEdges(space) << std::endl;
     std::cout << "connected" << lemon::connected(space) << std::endl;

     std::cout << "drawing" << std::endl;
     cliques::draw_graph canvas(space);
     cliques::make_energy_edge_colour_map A(stabilities);
     canvas.add_node_map(cliques::make_energy_edge_colour_map(stabilities));
     canvas.draw("spaces2.png");*/

    return 0;
}
