/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <stdlib.h>
#include <random>
#include <functional>
#include <list>

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

#include <cliques/algorithms/kernighan_lin.h>

namespace cliques {

// TODO - check num neighbours of space corresponds to algorithm exactly

/**
 @brief  Find all (single node moveset) neighbours of a partition
 */
template<typename G, typename P, typename NO>
bool will_move_break_partition(G &graph, const P &partition_const,
        NO &node_to_move) {
    typedef typename G::Node Node;
    typedef typename G::IncEdgeIt IncEdgeIt;

    P partition = partition_const;
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

 //TODO: This can return (and does return the same partition as a neighbour
 * Needs to be fixed.  This happens for example in a set of singletons, because
 * isolating a node would not break a partition
 */
template<typename G, typename P>
std::unordered_set<P, cliques::partition_hash, cliques::partition_equal> find_neighbours(
        G &graph,
        P const &partition) {
    typedef typename G::EdgeIt EdgeIt;
    typedef typename G::Node Node;

    std::unordered_set<P, cliques::partition_hash, cliques::partition_equal>
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
        if (will_move_break_partition(graph, partition, n1) == false) {
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

        if (will_move_break_partition(graph, partition, n2) == false) {
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

template<typename G, typename P, typename E, typename R>
P find_random_connected_neighbour(G &graph, P const &partition, E &edge,
        R &rand_engine) {
    typedef typename G::Node Node;

    Node n1 = graph.u(edge);
    Node n2 = graph.v(edge);
    int n1_id = graph.id(n1);
    int n2_id = graph.id(n2);
    int set_of_n1 = partition.find_set(n1_id);
    int set_of_n2 = partition.find_set(n2_id);
    bool are_in_same_set = (set_of_n1 == set_of_n2);
    std::unordered_set<P, cliques::partition_hash, cliques::partition_equal>
            neighbour_partitions;

    // Add partition with n1 isolated and in n2's set
    // Avoid using too much memory, destroy temp_partiton after use
    if (will_move_break_partition(graph, partition, n1) == false) {
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

    if (will_move_break_partition(graph, partition, n2) == false) {
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

    std::uniform_int_distribution<int> distribution(0,
            neighbour_partitions.size() - 1);

    int rand_neigh = distribution(rand_engine);
    auto set_itr = neighbour_partitions.begin();
    if (neighbour_partitions.size() > 0) {
        for (int i = 0; i < rand_neigh; ++i) {
            ++set_itr;
        }
        return *set_itr;
    } else {
        return partition;
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
        G &graph,
        std::unordered_set<P, cliques::partition_hash,
                cliques::partition_equal> &all_partitions, G &space) {

    typedef typename G::Node Node;
    typedef typename G::Edge Edge;
    typedef typename std::unordered_set<P, cliques::partition_hash,
            cliques::partition_equal> partition_set;
    typedef typename std::unordered_set<P, cliques::partition_hash,
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
        partition_set neighbour_partitions = cliques::find_neighbours(graph, *itr);

        for (partition_set_itr neigh_itr = neighbour_partitions.begin(); neigh_itr
                != neighbour_partitions.end(); ++neigh_itr) {

            // TODO: Hack! find_neighbours can return disconnected partition
            // Because of isolation
            // This discards if not in the set of all connected partitions
            if (all_partitions.find(*neigh_itr) == all_partitions.end()) {
                continue;
            }
            Node neighbour_node = partition_tofrom_Node.left.at(*neigh_itr);

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

template<typename G, typename DG, typename M>
void graph_to_transition_digraph(G &graph, std::vector<double> qualities,
        DG &transition_graph, M &transition_weights) {
    typedef typename G::NodeIt NodeIt;
    typedef typename G::IncEdgeIt IncEdgeIt;
    typedef typename G::Node Node;
    typedef typename DG::Arc Arc;

    int num_nodes = lemon::countNodes(graph);
    for (int i = 0; i < num_nodes; ++i) {
        transition_graph.addNode();
    }
    for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
        int base_id = graph.id(n1);

        double current_quality = qualities[base_id];
        double total_weight = 0.0;
        std::vector<Arc> unormalised_arcs;
        for (IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
            Node running = graph.oppositeNode(n1, e);
            int running_id = graph.id(running);
            double running_quality = qualities[running_id];
            if (running_quality > current_quality) {
                double weight_difference = running_quality - current_quality;
                total_weight += weight_difference;
                Arc new_arc = transition_graph.addArc(
                        transition_graph.nodeFromId(base_id),
                        transition_graph.nodeFromId(running_id));
                unormalised_arcs.push_back(new_arc);
                transition_weights.set(new_arc, weight_difference);
            }
        }
        for (typename std::vector<Arc>::iterator arc = unormalised_arcs.begin(); arc
                != unormalised_arcs.end(); ++arc) {
            double normalised_transition_weight = transition_weights[*arc]
                    / total_weight;
            transition_weights.set(*arc, normalised_transition_weight);
        }
    }
}

/**
 @brief  From a graph and vector of qualities create a directed transition graph
 The new graph has the same number of nodes, but only directed edges from one node
 to a node of greater quality, with the edge weight the probability of going
 to that node.

 The probability is weighted, as the difference between the qualities of the nodes
 / the total weight (positive with respect to the starting node)

 */
template<typename DG, typename M, typename NO>
void find_basins_depth_first(DG &transition_graph, M& transition_weights,
        double alpha, NO node, std::map<int, double> &basin_to_probabilities) {
    typedef typename DG::OutArcIt OutArcIt;
    typedef typename DG::InArcIt InArcIt;
    typedef typename DG::Node Node;

    for (InArcIt arc(transition_graph, node); arc != lemon::INVALID; ++arc) {
        Node lower_node = transition_graph.runningNode(arc);
        int lower_node_id = transition_graph.id(lower_node);
        double weight = transition_weights[arc];
        double new_alpha = weight * alpha;
        basin_to_probabilities[lower_node_id] += new_alpha;
        find_basins_depth_first(transition_graph, transition_weights,
                new_alpha, lower_node, basin_to_probabilities);
    }

}

template<typename G>
std::map<int, std::map<int, double> > compute_probabalistic_basins(G &graph,
        std::vector<double> qualities) {

    typedef typename lemon::SmartDigraph::NodeIt NodeIt;
    typedef typename lemon::SmartDigraph::OutArcIt OutArcIt;

    lemon::SmartDigraph transition_graph;
    lemon::SmartDigraph::ArcMap<double> transition_weights(transition_graph);
    graph_to_transition_digraph(graph, qualities, transition_graph,
            transition_weights);

    std::map<int, std::map<int, double> > all_basins;

    for (NodeIt node(transition_graph); node != lemon::INVALID; ++node) {
        int num_better_nodes = 0;
        for (OutArcIt arc(transition_graph, node); arc != lemon::INVALID; ++arc) {
            num_better_nodes++;
        }

        // I.e. only if this node is a maxima
        if (num_better_nodes == 0) {
            cliques::output("maxima");
            int node_id = transition_graph.id(node);
            double alpha = 1.0;
            std::map<int, double> basin_to_probabilities;
            find_basins_depth_first(transition_graph, transition_weights,
                    alpha, node, basin_to_probabilities);
            all_basins[node_id] = basin_to_probabilities;
        }
    }

    return all_basins;

}

/**
 @brief  From a graph and vector of qualities create a directed transition graph
 The new graph has the same number of nodes, but only directed edges from one node
 to a node of greater quality, with the edge weight the probability of going
 to that node.

 The probability is weighted, as the difference between the qualities of the nodes
 / the total weight (positive with respect to the starting node)

 */
template<typename G>
std::map<int, std::map<int, double> > compute_probabalistic_basins_new(
        G &graph, std::vector<double> qualities) {

    typedef typename lemon::SmartDigraph::NodeIt NodeIt;
    typedef typename lemon::SmartDigraph::Node Node;
    typedef typename lemon::SmartDigraph::OutArcIt OutArcIt;

    lemon::SmartDigraph transition_graph;
    lemon::SmartDigraph::ArcMap<double> transition_weights(transition_graph);
    graph_to_transition_digraph(graph, qualities, transition_graph,
            transition_weights);

    std::map<int, std::map<int, double> > all_basins;
    std::set<int> maxima;
    std::set<Node> maxima_nodes;
    lemon::SmartDigraph::NodeMap<int> real_outlinks_queue(transition_graph, 0);

    //find all maxima and insert into set
    for (NodeIt node(transition_graph); node != lemon::INVALID; ++node) {
        int num_better_nodes = 0;
        for (OutArcIt arc(transition_graph, node); arc != lemon::INVALID; ++arc) {
            num_better_nodes++;
        }
        real_outlinks_queue[node] = num_better_nodes;
        // I.e. only if this node is a maxima
        if (num_better_nodes == 0) {
            cliques::output("maxima");
            int node_id = transition_graph.id(node);
            maxima.insert(node_id);
            maxima_nodes.insert(node);
        }

    }

    for (std::set<int>::iterator itr = maxima.begin(); itr != maxima.end(); ++itr) {
        // for each maximum do message passing
        lemon::SmartDigraph::NodeMap<double> node_values(transition_graph, 0);
        int max_node_id = *itr;
        Node max_node = transition_graph.nodeFromId(max_node_id);
        node_values[max_node] = 1;
        std::set<Node> candidates(maxima_nodes);
        lemon::SmartDigraph::ArcMap<double> weights_to_max(transition_graph);
        lemon::mapCopy(transition_graph, transition_weights, weights_to_max);

        lemon::SmartDigraph::NodeMap<int> outlinks_queue(transition_graph);
        lemon::mapCopy(transition_graph, real_outlinks_queue, outlinks_queue);

        // as long as there are "unconverged" Nodes
        while (!candidates.empty()) {

            //loop over Nodes to be considered
            for (std::set<Node>::iterator candidate = candidates.begin(); candidate
                    != candidates.end(); ++candidate) {
                Node candidate_node = *candidate;

                // sum up outgoing link weight..
                for (lemon::SmartDigraph::OutArcIt outgoing_link(
                        transition_graph, candidate_node); outgoing_link
                        != lemon::INVALID; ++outgoing_link) {
                    node_values[candidate_node]
                            += weights_to_max[outgoing_link];
                }

                // multiply node value "down"..
                for (lemon::SmartDigraph::InArcIt incoming_link(
                        transition_graph, candidate_node); incoming_link
                        != lemon::INVALID; ++incoming_link) {

                    weights_to_max[incoming_link]
                            *= node_values[candidate_node];
                    Node downstream_node = transition_graph.source(
                            incoming_link);
                    outlinks_queue[downstream_node] -= 1;

                    // add potential nodes to candidate list
                    if (outlinks_queue[downstream_node] == 0) {
                        candidates.insert(downstream_node);
                    }

                }
                // remove node that passed on message
                candidates.erase(candidate_node);
            }
        }

        std::map<int, double> basin_to_probabilities;
        for (NodeIt node(transition_graph); node != lemon::INVALID; ++node) {
            int node_id = transition_graph.id(node);
            if (node_values[node]) {
                basin_to_probabilities[node_id] = node_values[node];
            }
        }

        all_basins[max_node_id] = basin_to_probabilities;
    }
    return all_basins;
}

/**
 @brief  Compute the probabalistic basins defined by kernginal lin
 algorithm
 */
template<typename G, typename M, typename QF, typename QFDIFF, typename MAP,
        typename COMM>
std::map<int, std::map<int, double> > compute_kernighan_lin_basins(G &graph,
        M &weights, QF compute_quality, QFDIFF compute_quality_diff, MAP &map,
        G &space, double time, COMM &all_partitions) {

    std::map<int, std::map<int, int> > node_to_basin_to_count;
    typedef cliques::VectorPartition VecPartition;
    const int TOTAL_COUNT_INDEX = -1; // Magic number to store total count
    int num_nodes = lemon::countNodes(graph);

    for (auto start_partition = all_partitions.begin(); start_partition
            != all_partitions.end(); ++start_partition) {

        std::vector<VecPartition> path_partitions;
        VecPartition output_partition(num_nodes);

        cliques::refine_partition_kernighan_lin_hijack(graph, weights,
                compute_quality, compute_quality_diff, path_partitions,
                *start_partition, output_partition);

        output_partition.normalise_ids();

        //Convert path_partitions into path_nodes
        std::vector<int> path_nodes;
        //        cliques::output("paths");
        for (auto partition = path_partitions.begin(); partition
                != path_partitions.end(); ++partition) {
            partition->normalise_ids();
            //            cliques::print_partition_line(*partition);
            auto p = map.left.find(*partition);
            if (p != map.left.end()) {
                auto node_id = space.id(p->second);
                path_nodes.push_back(node_id);
            }
        }

        //        cliques::output("start then output:");
        //        cliques::print_partition_line(*start_partition);
        //        cliques::print_partition_line(output_partition);

        auto p = map.left.find(output_partition);
        if (p != map.left.end()) {
            int optimum = space.id(p->second);
            path_nodes.push_back(space.id(map.left.at(*start_partition)));

            for (auto node = path_nodes.begin(); node != path_nodes.end(); ++node) {
                node_to_basin_to_count[*node][optimum]++;
                node_to_basin_to_count[*node][TOTAL_COUNT_INDEX]++;

                //                cliques::output(*node,optimum,node_to_basin_to_count[*node][optimum],node_to_basin_to_count[*node][TOTAL_COUNT_INDEX]);
            }
        } else {
            cliques::output("maxima is disconnected partition");
        }

    }

    std::map<int, std::map<int, double> > basin;

    // Convert into basin format
    for (auto node = node_to_basin_to_count.begin(); node
            != node_to_basin_to_count.end(); ++node) {
        int node_id = node->first;
        int total_count = node->second[TOTAL_COUNT_INDEX];

        for (auto basin_itr = node->second.begin(); basin_itr
                != node->second.end(); ++basin_itr) {
            int basin_id = basin_itr->first;
            if (basin_id != -1) {
                int basin_count = basin_itr->second;
                basin[basin_id][node_id] = basin_count / float(total_count);
            }
        }
    }

    return basin;
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
                        cliques::output(comm1_size, comm2_size, size_difference);
                    }
                }
            }
        }
    }
}

}
