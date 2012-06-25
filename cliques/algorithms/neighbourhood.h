#pragma once

namespace clq {

/**
 @brief  Find all (single node moveset) neighbours of a partition

 Finds neighbours of a partition where a neighbour is a partition
 which can be created by moving one node into an adjacent group
 or by isolating it into its own group.

 Basic algorithm: Iterate through edges, for each node u, v of edge:
 If moving the node (or isolation) would not break the partition:
 1. Isolate it
 2. If u and v are not in the same set, move u to v's set

 @param[in] all_partitons reference to unordered set of partitions
 @param[out] space output graph representing the space

 //TODO: This can return (and does) return the same partition as a neighbour
 * Needs to be fixed.  This happens for example in a set of singletons, because
 * isolating a node would not break a partition
 */
template<typename G, typename P>
std::unordered_set<P, clq::partition_hash, clq::partition_equal>
find_neighbours(
        G const &graph,
        P const &partition,
        bool allow_disconnected = true) {
    typedef typename G::EdgeIt EdgeIt;
    typedef typename G::Node Node;

    std::unordered_set<P, clq::partition_hash, clq::partition_equal>
            neighbour_partitions;

    for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
        Node n1 = graph.u(edge);
        Node n2 = graph.v(edge);
        int n1_id = graph.id(n1);
        int n2_id = graph.id(n2);
        int set_of_n1 = partition.find_set(n1_id);
        int set_of_n2 = partition.find_set(n2_id);
        bool are_in_same_set = (set_of_n1 == set_of_n2);

        // Add partition with n1 isolated and in n2's set
        // Avoid using too much memory, destroy temp_partiton after use
        if (allow_disconnected || will_move_break_partition(graph, partition, n1) == false) {
            P temp_partition = partition;
            temp_partition.isolate_node(n1_id);
            temp_partition.normalise_ids();
            neighbour_partitions.insert(temp_partition);

            if (!are_in_same_set) {
                P temp_partition = partition;
                temp_partition.add_node_to_set(n1_id, set_of_n2);
                temp_partition.normalise_ids();
                neighbour_partitions.insert(temp_partition);
            }
        }

        if (allow_disconnected || will_move_break_partition(graph, partition, n2) == false) {
            P temp_partition = partition;
            temp_partition.isolate_node(n2_id);
            temp_partition.normalise_ids();
            neighbour_partitions.insert(temp_partition);

            if (!are_in_same_set) {
                P temp_partition = partition;
                temp_partition.add_node_to_set(n2_id, set_of_n1);
                temp_partition.normalise_ids();
                neighbour_partitions.insert(temp_partition);
            }
        }
    }
    return neighbour_partitions;
}

/**
 @brief  Same as find_neighbours, but converts to stl vector
 */
template<typename G, typename P>
std::vector<P> find_neighbours_vec(
        G const &graph,
        P const &partition) {
    typedef typename G::EdgeIt EdgeIt;
    typedef typename G::Node Node;

    auto neighbour_partitions = find_neighbours(graph, partition);

    // Convert to vector
    std::vector<P> neighbour_partitions_vec;
    for (P neighbour : neighbour_partitions) {
        neighbour_partitions_vec.push_back(neighbour);
    }

    return neighbour_partitions_vec;
}

}