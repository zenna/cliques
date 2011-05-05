/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_SPACE_H
#define CLIQUES_SPACE_H

#include <boost/unordered_map.hpp>

#include <boost/bimap.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/multiset_of.hpp>

#include <cliques/graphhelpers.h>

namespace cliques {

/**
 @brief  Find all (single node moveset) neighbours of a partition

 Finds neighbours of a partition where a neighbour is a partition
 which can be created by moving one node into an adjacent group
 or by isolating it into its own group.

 @param[in] all_partitons reference to unordered set of partitions
 @param[out] space output graph representing the space
 */
template<typename G, typename P>
void find_neighbours(G &graph, P const &partition,
        boost::unordered_set<P, cliques::partition_hash,
                cliques::partition_equal> &neighbour_partitions) {
    P temp_partition = partition;

    typedef typename G::EdgeIt EdgeIt;

    for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
        int n1 = graph.id(graph.u(edge));
        int n2 = graph.id(graph.v(edge));
        int set_of_n1 = partition.find_set(n1);
        int set_of_n2 = partition.find_set(n2);

        // Add partition with n1 isolated and in n2's set
        // Avoid using too much memory, destroy temp_partiton after use
        {
            P temp_partition = partition;
            temp_partition.isolate_node(n1);
            temp_partition.normalise_ids();
            neighbour_partitions.insert(temp_partition);
        }
        {
            P temp_partition = partition;
            temp_partition.add_node_to_set(n1, set_of_n2);
            temp_partition.normalise_ids();
            neighbour_partitions.insert(temp_partition);
        }
        {
            P temp_partition = partition;
            temp_partition.isolate_node(n2);
            temp_partition.normalise_ids();
            neighbour_partitions.insert(temp_partition);
        }
        {
            P temp_partition = partition;
            temp_partition.add_node_to_set(n2, set_of_n1);
            temp_partition.normalise_ids();
            neighbour_partitions.insert(temp_partition);
        }
    }
}

/**
 @brief  Convert a set of partitions into a graph

 Given a set of partitions (perhaps found with cliques::find_connected_partitions),
 this function creates a graph with nodes representing partitions and edges created
 according to a moveset rule.  The default rule is the single node moveset
 see cliques::find_neighbours

 @param[in] all_partitons reference to unordered set of partitions
 @param[out] space output graph representing the space
 */
template<typename G, typename P>
boost::bimap<boost::bimaps::unordered_set_of<P, cliques::partition_hash,
        cliques::partition_equal>, boost::bimaps::set_of<typename G::Node> > create_space(
        G &graph, boost::unordered_set<P, cliques::partition_hash,
                cliques::partition_equal> &all_partitions, G &space) {

    typedef typename G::Node Node;
    typedef typename G::Edge Edge;
    typedef typename boost::unordered_set<P, cliques::partition_hash,
            cliques::partition_equal> partition_set;
    typedef typename boost::unordered_set<P, cliques::partition_hash,
            cliques::partition_equal>::iterator partition_set_itr;
    typedef boost::bimap<boost::bimaps::unordered_set_of<P,
            cliques::partition_hash, cliques::partition_equal>,
            boost::bimaps::set_of<Node> > Bimap;

    typedef typename Bimap::value_type bimap_value;
    Bimap partition_tofrom_Node;

    // Create node in space for each partition
    // And add to bimap relating node in space_map to partition
    for (partition_set_itr itr = all_partitions.begin(); itr
            != all_partitions.end(); ++itr) {
        Node temp_node = space.addNode();
        partition_tofrom_Node.insert(bimap_value(*itr, temp_node));
    }

    // For each partition, finds its neighbours and create edges
    // in space map
    for (partition_set_itr itr = all_partitions.begin(); itr
            != all_partitions.end(); ++itr) {

        Node current_node = partition_tofrom_Node.left.at(*itr);
        if (space.id(current_node) % 1000 == 0) {
            std::cout << space.id(current_node) << "\n";
        }
        partition_set neighbour_partitions;
        cliques::find_neighbours(graph, *itr, neighbour_partitions);

        for (partition_set_itr neigh_itr = neighbour_partitions.begin(); neigh_itr
                != neighbour_partitions.end(); ++neigh_itr) {

            // TODO: Hack! find_neighbours can return disconnected partition
            // This discards if not in the set of all connected partitions
            if (all_partitions.find(*neigh_itr) == all_partitions.end()) {
                continue;
            }
            Node neighbour_node = partition_tofrom_Node.left.at(*neigh_itr);

            if (current_node != neighbour_node) {
                Edge e = lemon::findEdge(space, current_node, neighbour_node);
                if (e == lemon::INVALID) {
                    Edge e = space.addEdge(current_node, neighbour_node);
                }
            }
        }
    }
    return partition_tofrom_Node;
}

template<typename G>
void create_disconnectivity_graph(G graph_landscape, G graph_dg) {
}

}

#endif
