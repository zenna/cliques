/* Copyright (c) Zenna Tavares, Michael Schaub - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <vector>
#include <map>
#include <limits>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#include <cliques/helpers/graphhelpers.h>
#include <cliques/quality_functions/internals/internals.h>
#include <cliques/quality_functions/internals/generators.h>


namespace cliques {
/**
 @brief  Refines a given partition using Kernighan Lin style movements

 This is a local algorithm which is based on modularity maximation outlined
 in Networks An Introduction, Newman page 373.  This in turn is based on
 an algorithm for finding graph bisection by Kernighan and Lin.

 It must take an initial partition, which may be random, but is typically the
 output of another optimisation such as louvain, upon which local optimisations
 are made.

 It works by looking at all pairs of nodes, and finding how the transition of node1
 into the community of n2 will change the quality (e.g. stability).  Over all
 pairs it executes the single best transition.  This process is repeated
 N times, while not moving nodes that have already been moved.
 Over the entire history of N transitions, it then finds the highest quality
 partition.  The algorithm is then recursed with the partition acting as the
 original partition as long as the quality continues to improve.

 Unlike louvain, this algorithm can make transitions which result in losses
 in quality in short run.

 @param[in]  graph     graph to find partition of
 @param[in]  compute_quality     partition quality function object
 @param[out]  std::vector<partitions>     optimal partitions, last in vector is overall best
 */
template<typename P, typename T, typename W, typename QF, typename QFDIFF>
double refine_partition_kernighan_lin(T &graph, W &weights, QF compute_quality,
        QFDIFF compute_quality_diff, P const &input_partition,
        P &output_partition) {
    typedef typename T::Node Node;
    typedef typename T::NodeIt NodeIt;
    typedef typename T::IncEdgeIt IncEdgeIt;

    // Initialisations..
    int num_nodes = lemon::countNodes(graph);
    P partition = input_partition;
    P buffer_partition = partition;
    auto internals = cliques::gen(compute_quality, graph, weights, partition);
    double minimum_improve = 0.000000001;
    double original_quality, best_quality, current_quality;
    original_quality = best_quality = current_quality = compute_quality(
            internals);

//    cliques::output("beginning", best_quality);
//    cliques::print_partition_line(partition);

    // keep track of nodes that have moved before
    std::set<Node> moved_nodes; //TODO change to unordered_set

    // Loop over number of nodes -- each node has to move once
    for (int i = 0; i < num_nodes; ++i) {
//        std::cout << "\nROUND  " << i << std::endl;
        // Initializations
        double best_gain = -std::numeric_limits<float>::max();
        Node node_to_move;
        unsigned int comm_to_move_to =-1;
        double absolute_gain;

        //		bool is_trapped_node = true;
        for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {

            // Don't move already moved nodes
            if (moved_nodes.find(n1) != moved_nodes.end()) {
                continue;
            }

            // get id of considered node and node community
            int n1_id = graph.id(n1);
            int n1_comm_id = partition.find_set(n1_id);

            //TODO check if there is a faster way of computing this as there are lots of updates involved here..
            // temporarily isolate node to compute the absolute gain via gain - loss

            isolate_and_update_internals(graph, weights, n1, internals,
                    partition);
            
            double isolation_loss = compute_quality_diff(internals, n1_comm_id,
                    n1_id);
            
            absolute_gain = -isolation_loss;
            if (absolute_gain > best_gain) {
                best_gain = absolute_gain;
                node_to_move = n1;
                int number_nodes = partition.element_count();
                std::vector<int> comm_existance(number_nodes, 0);
                for (int i = 0; i < partition.element_count(); ++i) {
                    int comm_id = partition.find_set(i);
                    if (comm_id != -1) {
                        comm_existance[comm_id]++;
                    }
                }
                // node was actually a singleton => gain = 0
                if(comm_existance[n1_comm_id] == 0){
                    best_gain =0;
                }
                for (int i = 0; i < partition.element_count(); ++i) {
                    if (comm_existance[i] == 0) {
                        comm_to_move_to = i;
                        break;
                    }
                }
            }

            // consider all possible other nodes
            // TODO: Don't we just have to look at all possible other communities the node could move to (at most all neighbouring nodes?)
            for (IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
                Node n2 = graph.oppositeNode(n1, e);
                int n2_id = graph.id(n2);
                int n2_comm_id = partition.find_set(n2_id);

                if (n1_comm_id == n2_comm_id) {
                    continue;
                }

                double gain =
                        compute_quality_diff(internals, n2_comm_id, n1_id);

                absolute_gain = gain - isolation_loss;

                // keep track of best possible move
                if (absolute_gain > best_gain) {
                    best_gain = absolute_gain;
                    node_to_move = n1;
                    comm_to_move_to = n2_comm_id;
                }
            }
            // Put node back - we don't want to move yet
            insert_and_update_internals(graph, weights, n1, internals,
                    partition, n1_comm_id);
        }

        // TODO: check if this can be done more efficient, see above
        // move node from old community to other
        isolate_and_update_internals(graph, weights, node_to_move, internals,
                partition);
        insert_and_update_internals(graph, weights, node_to_move, internals,
                partition, comm_to_move_to);
        // keep track of moved nodes
        moved_nodes.insert(node_to_move);

        // keep track of quality
        current_quality = current_quality + best_gain;


        if (current_quality > best_quality) {
//            cliques::output(current_quality, best_quality);
            buffer_partition = partition;
            best_quality = current_quality;
        }
    }



    double total_quality_improvement = best_quality - original_quality;
//     cliques::print_partition_line(buffer_partition);
//     cliques::output(current_quality,total_quality_improvement, best_quality, original_quality);

    // TODO maybe it is better to have an iterative instead of recursive implementation as the stack might grow large unnecessarily
    if (total_quality_improvement > minimum_improve) {
//        cliques::output("end",best_quality);
//        cliques::print_partition_line(buffer_partition);
        return total_quality_improvement + refine_partition_kernighan_lin(
                graph, weights, compute_quality, compute_quality_diff,
                buffer_partition, output_partition);
    } else {
        output_partition = buffer_partition;
        return total_quality_improvement;
    }
}

// ------------------------------------------------------------------------

/**
 @brief  Modified Kernighan Lin for experimentation
 // return a list of partitions
 */
template<typename P, typename T, typename W, typename QF, typename QFDIFF>
double refine_partition_kernighan_lin_hijack(T &graph, W &weights, QF compute_quality,
        QFDIFF compute_quality_diff, std::vector<P> &path_partitions, P const &input_partition,
        P &output_partition) {
    typedef typename T::Node Node;
    typedef typename T::NodeIt NodeIt;
    typedef typename T::IncEdgeIt IncEdgeIt;

    // Initialisations..
    int num_nodes = lemon::countNodes(graph);
    P partition = input_partition;
    P buffer_partition = partition;
    auto internals = cliques::gen(compute_quality, graph, weights, partition, partition);
    double minimum_improve = 0.000000001;
    double original_quality, best_quality, current_quality;
    original_quality = best_quality = current_quality = compute_quality(
            internals);

//    cliques::output("beginning", best_quality);
//    cliques::print_partition_line(partition);

    // keep track of nodes that have moved before
    std::set<Node> moved_nodes; //TODO change to unordered_set

    // Loop over number of nodes -- each node has to move once
    for (int i = 0; i < num_nodes; ++i) {
        double best_gain = -std::numeric_limits<float>::max();
        Node node_to_move;
        unsigned int comm_to_move_to = -1;
        double absolute_gain;

        //      bool is_trapped_node = true;
        for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {

            // Don't move already moved nodes
            if (moved_nodes.find(n1) != moved_nodes.end()) {
                continue;
            }

            // get id of considered node and node community
            int n1_id = graph.id(n1);
            int n1_comm_id = partition.find_set(n1_id);

            //TODO check if there is a faster way of computing this as there are lots of updates involved here..
            // temporarily isolate node to compute the absolute gain via gain - loss

            isolate_and_update_internals(graph, weights, n1, internals,
                    partition);

            double isolation_loss = compute_quality_diff(internals, n1_comm_id,
                    n1_id);

            absolute_gain = -isolation_loss;
            if (absolute_gain > best_gain) {
                best_gain = absolute_gain;
                node_to_move = n1;
                int number_nodes = partition.element_count();
                std::vector<int> comm_existance(number_nodes, 0);
                for (int i = 0; i < partition.element_count(); ++i) {
                    int comm_id = partition.find_set(i);
                    if (comm_id != -1) {
                        comm_existance[comm_id]++;
                    }
                }
                // node was actually a singleton => gain = 0
                if(comm_existance[n1_comm_id] == 0){
                    best_gain =0;
                }
                for (int i = 0; i < partition.element_count(); ++i) {
                    if (comm_existance[i] == 0) {
                        comm_to_move_to = i;
                        break;
                    }
                }
            }

            // consider all possible other nodes
            // TODO: Don't we just have to look at all possible other communities the node could move to (at most all neighbouring nodes?)
            for (IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
                Node n2 = graph.oppositeNode(n1, e);
                int n2_id = graph.id(n2);
                int n2_comm_id = partition.find_set(n2_id);

                if (n1_comm_id == n2_comm_id) {
                    continue;
                }

                double gain =
                        compute_quality_diff(internals, n2_comm_id, n1_id);

                absolute_gain = gain - isolation_loss;

                // keep track of best possible move
                if (absolute_gain > best_gain) {
                    best_gain = absolute_gain;
                    node_to_move = n1;
                    comm_to_move_to = n2_comm_id;
                }
            }
            // Put node back - we don't want to move yet
            insert_and_update_internals(graph, weights, n1, internals,
                    partition, n1_comm_id);
        }

        // TODO: check if this can be done more efficient, see above
        // move node from old community to other
        isolate_and_update_internals(graph, weights, node_to_move, internals,
                partition);
        insert_and_update_internals(graph, weights, node_to_move, internals,
                partition, comm_to_move_to);

        path_partitions.push_back(partition);
        // keep track of moved nodes
        moved_nodes.insert(node_to_move);

        // keep track of quality
        current_quality = current_quality + best_gain;


        if (current_quality > best_quality) {
//            cliques::output(current_quality, best_quality);
            buffer_partition = partition;
            best_quality = current_quality;
        }
    }

    double total_quality_improvement = best_quality - original_quality;
//     cliques::print_partition_line(buffer_partition);
//     cliques::output(current_quality,total_quality_improvement, best_quality, original_quality);

    // TODO maybe it is better to have an iterative instead of recursive implementation as the stack might grow large unnecessarily
    if (total_quality_improvement > minimum_improve) {
//        cliques::output("end",best_quality);
//        cliques::print_partition_line(buffer_partition);
        return total_quality_improvement + refine_partition_kernighan_lin(
                graph, weights, compute_quality, compute_quality_diff,
                buffer_partition, output_partition);
    } else {
        output_partition = buffer_partition;
        return total_quality_improvement;
    }
}

}
