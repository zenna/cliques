/* Copyright (c) Zenna Tavares-zennatavares@gmail.com, Michael Schaub 2010-2012 */
#pragma once

#include <vector>
#include <limits>
#include <map>
#include <cliques/structures/vector_partition.h>

namespace clq {

template<typename G>
clq::VectorPartition community_to_partition(G &graph,
        std::vector<int> comm, int outside_comm_set_id) {
    int num_nodes = lemon::countNodes(graph);
    clq::VectorPartition partition(num_nodes, outside_comm_set_id);
    for (auto node = comm.begin(); node != comm.end(); ++node) {
        partition.add_node_to_set(*node, 1);
    }
    //    partition.normalise_ids();
    return partition;
}

/**
 @brief  Will removal of node_to_remove from graph cause community comm to be disconnected?
 */
template<typename G, typename NO>
bool will_removal_break_community(G &graph, std::vector<int> comm,
        NO &node_to_remove) {
    clq::VectorPartition partition = community_to_partition(graph, comm, 0);
    return clq::will_move_break_partition(graph, partition, node_to_remove);
}

template<typename G>
std::vector<std::vector<int> > find_community_neighbours(G &graph,
        std::vector<int> comm) {

    std::set<int> forbidden;
    std::vector<std::vector<int> > neighbours;

    // Handle special case of single node community
    if (comm.size() == 1) {
        forbidden.insert(comm[0]);
    }
    for (unsigned int i = 0; i < comm.size(); ++i) { // Better if I could just iterate thro
        int node_id = comm[i];
        auto node = graph.nodeFromId(node_id);

        // If removal won't split community
        if (forbidden.find(node_id) == forbidden.end()
                && will_removal_break_community(graph, comm, node) == false) {
            std::vector<int> contracted_comm = comm;
            auto location = std::find(contracted_comm.begin(),
                    contracted_comm.end(), node_id);
            contracted_comm.erase(location);
            forbidden.insert(node_id);
            std::sort(contracted_comm.begin(), contracted_comm.end());
            //		    clq::output("WHAAC");
            //		    clq::print_collection(contracted_comm);
            neighbours.push_back(contracted_comm);
        }

        for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
            auto neigh_node = graph.oppositeNode(node, e);
            int neigh_node_id = graph.id(neigh_node);

            // If node is on periphery
            if (std::find(comm.begin(), comm.end(), neigh_node_id)
                    == comm.end() && forbidden.find(neigh_node_id)
                    == forbidden.end()) {
                forbidden.insert(neigh_node_id);
                std::vector<int> new_comm = comm;
                new_comm.push_back(neigh_node_id);
                std::sort(new_comm.begin(), new_comm.end());

                neighbours.push_back(new_comm);//neighbours.push_back()
            }
        }
    }
    return neighbours;
}

/**
 @brief  Find maxima

 */
template<typename G, typename M, typename QF>
bool is_community_maxima(G &graph, M &weights, std::vector<int> comm,
        QF compute_quality, double &stability) {
    clq::VectorPartition comm_partition = community_to_partition(graph,
            comm, 0);
    auto internals = clq::gen(compute_quality, graph, weights,
            comm_partition);
    double current_quality = compute_quality(internals, 1, 2.3);
    double best_quality = -std::numeric_limits<float>::max();
    stability = current_quality;

    auto neighs = find_community_neighbours(graph, comm);
    for (auto itr = neighs.begin(); itr != neighs.end(); ++itr) {
        clq::VectorPartition neigh_partition = community_to_partition(
                graph, *itr, 0);
        auto neigh_internals = clq::gen(compute_quality, graph, weights,
                neigh_partition);
        double neigh_quality = compute_quality(neigh_internals, 1, 2.3);

        if (neigh_quality > best_quality) {
            best_quality = neigh_quality;
        }
    }

    // Neutrality?
    if (current_quality > best_quality) {
        return true;
    } else {
        return false;
    }

}

}
