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

#include <cliques/algorithms/space.h>
#include <cliques/structures/vector_partition.h>

int main() {
    lemon::SmartGraph orange_graph;
    lemon::SmartGraph::EdgeMap<int> weights(orange_graph);

    typedef cliques::VectorPartition VecPartition;
    cliques::read_edgelist_weighted(
            "/home/zenna/repos/graph-codes/cliques/data/graphs/renaud_n16.edj",
            orange_graph, weights);

    /*std::cout << "num_nodes" << lemon::countNodes(orange_graph) << std::endl;
     std::cout << "num_edges" << lemon::countEdges(orange_graph) << std::endl;
     std::cout << "connected" << lemon::connected(orange_graph) << std::endl;*/

    //std::cout << "Finding connected partitions" << std::endl;
    boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> all_partitions;
    cliques::find_connected_partitions(orange_graph, all_partitions);

    std::vector<double> markov_times = { 1.0 };
    double current_markov_time = 1.0;
    cliques::find_weighted_linearised_stability func(markov_times);
    int num_partitions = 0;

    double stability = 0.0;
    /*for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
        VecPartition p = *itr;
        cliques::Internals internals(orange_graph, weights, p);
        stability = func(internals);
        std::cout << stability << std::endl;
    }*/

    std::cout << "Computing Space" << std::endl;
    lemon::SmartGraph space;
    auto map = cliques::create_space(orange_graph, all_partitions, space);

    std::vector<VecPartition> optimal_partitions;
    stability
            = cliques::find_optimal_partition_louvain_with_gain<VecPartition>(
                    space, weights,
                    cliques::find_weighted_linearised_stability(markov_times),
                    cliques::linearised_stability_gain_louvain(
                            current_markov_time), optimal_partitions);

    VecPartition best_partition = optimal_partitions.back();
    best_partition.normalise_ids();
    std::cout << "num_sets: " << *std::max_element(
            best_partition.partition_vector.begin(),
            best_partition.partition_vector.end());

    /*std::vector<float> stabilities;
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
