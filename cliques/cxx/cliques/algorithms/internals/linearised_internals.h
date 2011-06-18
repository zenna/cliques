#pragma once

#include <cliques/algorithms/stability.h>

namespace cliques {
/**
 @brief  isolate a node into its singleton set & update internals
 */
template<typename G, typename M, typename I, typename P>
void isolate_and_update_internals(G &graph, M &weights, typename G::Node node,
        I &internals, P &partition) {
    int node_id = graph.id(node);
    int comm_id = partition.find_set(node_id);

    // if node is already isolated
    if (internals.comm_w_in[comm_id] != 0) {
        double summed_weight_inside = 0;
        double summed_weight_outside = 0;
        for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
            double edge_weight = weights[e];
            typename G::Node opposite_node = graph.oppositeNode(node, e);
            int comm_node = partition.find_set(graph.id(opposite_node));
            if (comm_node == comm_id) {
                summed_weight_inside += edge_weight;
            } else {
                summed_weight_outside += edge_weight;
            }
        }

        // When isolate a node, its old set looses the edges that node had to its outside set
        // but gains new edges it makes to inside the set
        internals.comm_to_sum_inc_weight[comm_id]
                = internals.comm_to_sum_inc_weight[comm_id]
                        - summed_weight_outside + summed_weight_inside;
        //            internals.comm_to_sum_inc_weight[new_comm_id] = summed_weight_outside + summed_weight_inside;
        //            internals.comm_w_tot[new_comm_id] = summed_weight_outside + summed_weight_inside;
        //            internals.comm_w_in[new_comm_id] = 0;

        internals.comm_w_tot[comm_id] -= internals.node_to_w[node_id];
        internals.comm_w_in[comm_id] -= 2 * summed_weight_inside
                + find_weight_selfloops(graph, weights, node);
        partition.isolate_node(node_id);

    }
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename G, typename M, typename I, typename P>
void insert_and_update_internals(G &graph, M &weights, typename G::Node node,
        I &internals, P &partition, int best_comm) {
    int node_id = graph.id(node);

    double summed_weight_inside = 0;
    double summed_weight_outside = 0;
    for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
        double edge_weight = weights[e];
        typename G::Node opposite_node = graph.oppositeNode(node, e);
        int comm_node = partition.find_set(graph.id(opposite_node));
        if (comm_node == best_comm) {
            summed_weight_inside += edge_weight;
        } else {
            summed_weight_outside += edge_weight;
        }
    }

    // Handle case if adding back to previous isolate comm_id
    // Don't want to accumulate comm_w_tot
    if (internals.comm_w_in[best_comm] == 0) {
        internals.comm_w_tot[best_comm] = internals.node_to_w[node_id];
    } else {
        internals.comm_w_tot[best_comm] += internals.node_to_w[node_id];
        internals.comm_w_in[best_comm] += 2
                * summed_weight_inside
                + find_weight_selfloops(graph, weights, node);
        internals.comm_to_sum_inc_weight[best_comm]
                = internals.comm_to_sum_inc_weight[best_comm]
                        - summed_weight_outside + summed_weight_inside;
    }
    partition.add_node_to_set(node_id, best_comm);
}

struct LinearisedInternals {
    typedef lemon::RangeMap<double> range_map;
    unsigned int num_nodes;
    double two_m;
    range_map node_to_w;
    range_map comm_w_tot;
    range_map comm_w_in;
    std::map<int, double> comm_to_sum_inc_weight;

    void print() {
        cliques::output("num_nodes:", num_nodes, "two_m:", two_m);
        cliques::output("node_to_w");
        for (int i = 0; i < node_to_w.size(); ++i) {
            std::cout << node_to_w[i] << " ";
        }
        std::cout << std::endl;
        cliques::output("comm_w_tot");
        for (int i = 0; i < comm_w_tot.size(); ++i) {
            std::cout << comm_w_tot[i] << " ";
        }
        std::cout << std::endl;
        cliques::output("comm_w_in");
        for (int i = 0; i < comm_w_in.size(); ++i) {
            std::cout << comm_w_in[i] << " ";
        }
        std::cout << std::endl;
        cliques::output("comm_to_sum_inc_weight");
        for (auto itr = comm_to_sum_inc_weight.begin(); itr
                != comm_to_sum_inc_weight.end(); ++itr) {
            std::cout << itr->first << " -> " << itr->second << std::endl;
        }
        std::cout << std::endl;
    }

    template<typename G, typename M>
    LinearisedInternals(G &graph, M &weights) :
        num_nodes(lemon::countNodes(graph)), node_to_w(num_nodes, 0),
                comm_w_tot(num_nodes, 0), comm_w_in(num_nodes, 0) {
        two_m = 2 * find_total_weight(graph, weights);
        for (unsigned int i = 0; i < num_nodes; ++i) {
            typename G::Node temp_node = graph.nodeFromId(i);
            comm_w_tot[i] = node_to_w[i] = find_weighted_degree(graph, weights,
                    temp_node);
            comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node);
        }
    }

    template<typename G, typename M, typename P>
    LinearisedInternals(G &graph, M &weights, P &partition) :
        num_nodes(lemon::countNodes(graph)), node_to_w(num_nodes, 0),
                comm_w_tot(num_nodes, 0), comm_w_in(num_nodes, 0) {
        two_m = 2 * find_total_weight(graph, weights);

        typedef typename G::EdgeIt EdgeIt;
        // find internal statistics based on graph, weights and partitions
        // consider all edges
        for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
            int node_u_id = graph.id(graph.u(edge));
            int node_v_id = graph.id(graph.v(edge));

            // this is to distinguish within community weight with total weight
            int comm_of_node_u = partition.find_set(node_u_id);
            int comm_of_node_v = partition.find_set(node_v_id);

            // weight of edge
            double weight = weights[edge];

            // add weight to node weight
            node_to_w[node_u_id] += weight;
            node_to_w[node_v_id] += weight;
            // add weight to total weight of community
            comm_w_tot[comm_of_node_u] += weight;
            comm_w_tot[comm_of_node_v] += weight;
            if (comm_of_node_u == comm_of_node_v) {
                // in case the weight stems from within the community add to internal weights
                comm_w_in[comm_of_node_u] += 2 * weight;
            } else {
                comm_to_sum_inc_weight[comm_of_node_u]
                        = comm_to_sum_inc_weight[comm_of_node_u] + 1;
                comm_to_sum_inc_weight[comm_of_node_v]
                        = comm_to_sum_inc_weight[comm_of_node_v] + 1;
            }
        }
    }
};

template<typename G, typename M>
cliques::LinearisedInternals gen(find_weighted_linearised_stability, G &graph,
        M &weights) {
    LinearisedInternals internals(graph, weights);
    return internals;
}

template<typename G, typename M, typename P>
cliques::LinearisedInternals gen(find_weighted_linearised_stability, G &graph,
        M &weights, P &partition) {
    LinearisedInternals internals(graph, weights, partition);
    return internals;
}
}
