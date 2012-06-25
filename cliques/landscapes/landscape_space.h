/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <stdlib.h>
#include <random>
#include <functional>
#include <list>

#include <unordered_set>
#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/unordered_set.hpp>
#include <boost/bimap/multiset_of.hpp>

#include <cliques/helpers/lists.h>
#include <cliques/helpers/graphhelpers.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/structures/common.h>
#include "cliques/algorithms/neighbourhood.h"

#include <cliques/algorithms/kernighan_lin.h>

namespace clq {

// TODO - check num neighbours of space corresponds to algorithm exactly

/**
 @brief  Will moving node_to_move break the partition?
 */
template<typename G, typename P, typename NO>
bool will_move_break_partition(G &graph, P partition,
        NO &node_to_move) {
    typedef typename G::Node Node;
    typedef typename G::IncEdgeIt IncEdgeIt;

    int node_to_move_id = graph.id(node_to_move);
    std::set<int> neighs_in_same_comm;
    int node_to_move_comm = partition.find_set(node_to_move_id);

    // Find neighbours of node in the same community
    for (IncEdgeIt e(graph, node_to_move); e != lemon::INVALID; ++e) {
        Node neigh_node = graph.runningNode(e);
        int neigh_node_id = graph.id(neigh_node);
        if (node_to_move_comm == partition.find_set(neigh_node_id)) {
            neighs_in_same_comm.insert(neigh_node_id);
        }
    }

    // If node is already isolated
    if (neighs_in_same_comm.size() == 0) {
        return false;
    }

    partition.isolate_node(node_to_move_id);
    std::set<int> seen_nodes;
    std::list<int> queue(1, *neighs_in_same_comm.begin());
    // Do breadth first search from a neighbour
    // Until either all neighs_in_same_comm found
    // or nowhere left to search
    while (queue.size() != 0) {
        int source_id = queue.back();
        Node source = graph.nodeFromId(source_id);
        queue.pop_back();
        // TODO: This only needs to be done once, for first source
        // the rest is unncessary computation
        seen_nodes.insert(source_id);

        // If you remove a node once it has been seen then this
        // Will always go to zero
        // If we don't remove it
        auto node_itr = neighs_in_same_comm.find(source_id);
        if (node_itr != neighs_in_same_comm.end()) {
            neighs_in_same_comm.erase(node_itr);
        }

        for (IncEdgeIt e(graph, source); e != lemon::INVALID; ++e) {
            Node neigh_node = graph.runningNode(e);
            int neigh_node_id = graph.id(neigh_node);

            // Only look at neighbours within same community
            if (partition.find_set(neigh_node_id) == node_to_move_comm) {
                if (seen_nodes.find(neigh_node_id) == seen_nodes.end()) {
                    queue.push_back(neigh_node_id);
                    seen_nodes.insert(neigh_node_id);
                } else {
                    continue;
                }

                auto neigh_itr = neighs_in_same_comm.find(neigh_node_id);
                if (neigh_itr != neighs_in_same_comm.end()) {
                    neighs_in_same_comm.erase(neigh_itr);
                }
            }
        }
        if (neighs_in_same_comm.size() == 0) {
            // Undo isolation
            partition.add_node_to_set(node_to_move_id, node_to_move_comm);
            return false;
        }
    }

    partition.add_node_to_set(node_to_move_id, node_to_move_comm);
    return true;
}

/**
 @brief  Convert a set of partitions into a graph

 Given a set of partitions (perhaps found with clq::find_connected_partitions),
 this function creates a graph with nodes representing partitions and edges created
 according to a moveset rule.  The default rule is the single node moveset
 see clq::find_neighbours

 @param[in] all_partitons reference to unordered set of partitions
 @param[out] space output graph representing the space
 */
template<typename G, typename P>
boost::bimap<boost::bimaps::unordered_set_of<P, clq::partition_hash,
        clq::partition_equal>, boost::bimaps::set_of<typename G::Node> >
create_space(
    G const &graph,
    std::unordered_set<P, clq::partition_hash,
        clq::partition_equal> const &all_partitions,
    G &space) {

    typedef typename G::Node Node;
    typedef typename G::Edge Edge;
    typedef typename std::unordered_set<P, clq::partition_hash,
            clq::partition_equal> partition_set;
    typedef boost::bimap<boost::bimaps::unordered_set_of<P,
            clq::partition_hash, clq::partition_equal>,
            boost::bimaps::set_of<Node> > Bimap;


    typedef typename Bimap::value_type bimap_value;
    Bimap partition_tofrom_Node;

    // Create node in space for each partition
    // And add to bimap relating node in space_map to partition
    for (P const &partition : all_partitions) {
        Node temp_node = space.addNode();
        partition_tofrom_Node.insert(bimap_value(partition, temp_node));
    }

    // Find neighbour of each partition and make edges in space map
    for (P const &partition : all_partitions) {

        Node current_node = partition_tofrom_Node.left.at(partition);
        if (space.id(current_node) % 1000 == 0) {
            std::cout << space.id(current_node) << "\n";
        }
        partition_set neighbour_partitions = clq::find_neighbours(graph, partition);

        for (P const &neigh_partition : neighbour_partitions) {

            // TODO: Hack! find_neighbours can return disconnected partition
            // This discards if not in the set of all connected partitions
            if (all_partitions.find(neigh_partition) == all_partitions.end()) {
                continue;
            }
            Node neighbour_node = partition_tofrom_Node.left.at(neigh_partition);

            if (current_node != neighbour_node) {
                Edge e = lemon::findEdge(space, current_node, neighbour_node);
                if (e == lemon::INVALID) {
                    space.addEdge(current_node, neighbour_node);
                }
            }
        }
    }
    return partition_tofrom_Node;
}

/**
 @brief  Convert a set of partitions into a graph

 Given a set of partitions (perhaps found with clq::find_connected_partitions),
 this function creates a graph with nodes representing partitions and edges created
 according to a moveset rule.  The default rule is the single node moveset
 see clq::find_neighbours

 @param[in] all_partitons reference to unordered set of partitions
 @param[out] space output graph representing the space
 */
template<typename G>
void create_hasse_community_space(G &graph,
        std::vector<std::vector<int> > &communities, G &space) {

    typedef typename G::Node Node;

    int num_nodes = communities.size();
    for (int i = 0; i < num_nodes; ++i) {
        space.addNode();
    }

    std::vector<int> buffer(20, 0);

    for (int i = 0; i < num_nodes - 1; ++i) {
        for (int j = i + 1; j < num_nodes; ++j) {
            auto comm1 = communities[i];
            auto comm2 = communities[j];
            int comm1_size = comm1.size();
            int comm2_size = comm2.size();
            int size_difference = comm1_size - comm2_size;

            // TODO: Does this really work here??
            // Communities of the same size cannot contain one another
            if (std::abs(size_difference) <= 13) {
                std::vector<int> &larger_comm = comm1_size > comm2_size ? comm1
                        : comm2;
                std::vector<int> &smaller_comm =
                        comm1_size > comm2_size ? comm2 : comm1;

                if (does_set_contain_other(larger_comm, smaller_comm, buffer)) {
                    Node u = graph.nodeFromId(i);
                    Node v = graph.nodeFromId(j);
                    graph.addEdge(u, v);
                    if (size_difference >= 11) {
                        clq::output(comm1_size, comm2_size, size_difference);
                    }
                }
            }
        }
    }
}

}
