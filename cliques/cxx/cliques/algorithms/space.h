/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <random>
#include <functional>
#include <list>

#include <boost/unordered_map.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/multiset_of.hpp>

#include <cliques/graphhelpers.h>
#include <cliques/structures/vector_partition.h>

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
void find_neighbours(
        G &graph,
        P const &partition,
        boost::unordered_set<P, cliques::partition_hash,
                cliques::partition_equal> &neighbour_partitions) {
    typedef typename G::EdgeIt EdgeIt;
    typedef typename G::Node Node;

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
    boost::unordered_set<P, cliques::partition_hash, cliques::partition_equal>
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
        boost::unordered_set<P, cliques::partition_hash,
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

/**
 @brief  Uniformly sample partition space using Metropolis-Hastings
 // TODO make type partition type independent
 */
//template<typename G, typename S, typename Logger>
//void sample_uniform_metropolis(G &graph, int num_samples,
//        int num_steps_per_sample, S &sampled_partitions, Logger &logger) {
//    typedef typename boost::unordered_set<cliques::VectorPartition,
//            cliques::partition_hash, cliques::partition_equal> partition_set;
//
//    int num_sampled = 0;
//    int num_steps = 0;
//    int num_nodes = lemon::countNodes(graph);
//    cliques::VectorPartition current_partition(num_nodes);
//    current_partition.initialise_as_singletons();
//
//    std::uniform_real_distribution<> real_distribution(0,1);
//    std::mt19937 m_engine;
//    std::mt19937 engine; // Mersenne twister MT19937
//
//    while (true) {
//        partition_set neigh_partitions;
//
//
//
//        std::vector<int>  debug_partition_vector = {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0};
//        cliques::VectorPartition debug_partition(debug_partition_vector);
//
//        if (current_partition == debug_partition) {
//            cliques::output("HOLD UP A MINUTE");
//        }
//
//        cliques::find_neighbours(graph, current_partition, neigh_partitions);
//        int num_current_neighs = neigh_partitions.size();
//
//        cliques::output("\nThe following partition", num_current_neighs);
//        cliques::print_partition_line(current_partition);
//        cliques::output("has", num_current_neighs, "neighbours");
//
//
//        for (auto a = neigh_partitions.begin(); a != neigh_partitions.end(); ++a) {
//            cliques::print_partition_line(*a);
//        }
//        std::uniform_int_distribution<int> distribution(0,
//                num_current_neighs-1);
//
//        while (true) {
//            int rand_neigh = distribution(engine);
//
//            auto set_itr = neigh_partitions.begin();
//            for (int i = 0; i < rand_neigh; ++i) {
//                ++set_itr;
//            }
//            cliques::VectorPartition proposed_partition = *set_itr;
//            logger.log(proposed_partition);
//            partition_set proposed_neighs;
//            //cliques::print_partition_list(proposed_partition);
//            cliques::find_neighbours(graph, proposed_partition, proposed_neighs);
//            int num_proposed_neighs = proposed_neighs.size();
//
//            //Metropolis-Hastings acceptance alpha
//            double alpha = real_distribution(m_engine);
//            double rand_real_num = double(num_current_neighs) / double(num_proposed_neighs);
//
//            if (alpha < rand_real_num) {
//                num_steps++;
//                current_partition = proposed_partition;
//                break;
//            }
//        }
//        if (num_steps % num_steps_per_sample == 0) {
//        	current_partition.normalise_ids();
//            sampled_partitions.insert(current_partition);
//            logger.log(current_partition);
//            num_sampled++;
//            //cliques::output(num_sampled, sampled_partitions.size());
//        }
//        if (num_sampled == num_samples) {
//            break;
//        }
//    }
//}

/**
 @brief  Uniformly sample partition space using Metropolis-Hastings
 // TODO make type partition type independent
 */
//template<typename G, typename P, typename Logger>
//void sample_uniform_biased(G &graph, int num_samples,
//        int num_steps_per_sample, boost::unordered_set<P, cliques::partition_hash,
//        cliques::partition_equal> &sampled_partitions, Logger &logger) {
//    typedef typename boost::unordered_set<P,
//            cliques::partition_hash, cliques::partition_equal> partition_set;
//    typedef typename G::Edge Edge;
//
//    int num_sampled = 0;
//    int num_steps = 0;
//    int num_nodes = lemon::countNodes(graph);
//    P current_partition(num_nodes);
//    current_partition.initialise_as_singletons();
//
//    std::mt19937 rand_engine; // Mersenne twister MT19937
//    std::vector<Edge> edges;
//
//    for (typename G::EdgeIt e(graph); e != lemon::INVALID; ++e) {
//    	edges.push_back(e);
//    }
//
//    std::uniform_int_distribution<int> distribution(0,
//            edges.size()-1);
//
//    while (true) {
//        Edge rand_edge = edges[distribution(rand_engine)];
//    	P new_partition = find_random_connected_neighbour(graph, current_partition, rand_edge, rand_engine);
//
//    	// Try again if no neighbour found
//    	if (new_partition == current_partition) {
//        	continue;
//        }
//        current_partition = new_partition;
//    	if (num_steps % num_steps_per_sample == 0) {
//        	current_partition.normalise_ids();
//            sampled_partitions.insert(current_partition);
//            logger.log(current_partition);
//            num_sampled++;
//        }
//        ++num_steps;
//        if (num_sampled == num_samples) {
//            break;
//        }
//    }
//}

// TODO
//Convert Graph To Transition Graph
// Find Basin Probabilities and write to file
// Is Partition Maxima
// Create communities space graph
// Create Hasse diagram type graph:
//
//template <typename G, typename DG>
//void convert_graph_to_transition_graph(G &graph, std::vector<dobule> energies, DG &transition_graph) {
//	lemon::INVALID;
//}
//
//void wh
//
//}

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
        for (auto arc = unormalised_arcs.begin(); arc != unormalised_arcs.end(); ++arc) {
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
void find_basins_depth_first(DG &transition_graph, M& transition_weights, double alpha,
        NO node, std::map<int, double> &basin_to_probabilities) {
    typedef typename DG::OutArcIt OutArcIt;
    typedef typename DG::InArcIt InArcIt;
    typedef typename DG::Node Node;

    for (InArcIt arc(transition_graph, node); arc != lemon::INVALID; ++arc) {
        Node lower_node = transition_graph.runningNode(arc);
        int lower_node_id = transition_graph.id(lower_node);
        double weight = transition_weights[arc];
        double new_alpha = weight * alpha;
        basin_to_probabilities[lower_node_id] += new_alpha;
        find_basins_depth_first(transition_graph, transition_weights, new_alpha, lower_node,
                basin_to_probabilities);
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

    std::map<int, std::map<int, double> >all_basins;

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
            find_basins_depth_first(transition_graph, transition_weights, alpha, node, basin_to_probabilities);
            all_basins[node_id] = basin_to_probabilities;
        }
    }

    return all_basins;

}

}
