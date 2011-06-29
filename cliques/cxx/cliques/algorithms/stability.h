/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/graphhelpers.h>
#include <cliques/structures/disjointset.h>

namespace cliques {
/**
 @brief  Functor for finding stability of weighted graph
 */
struct find_weighted_linearised_stability {
    double markov_time;

    find_weighted_linearised_stability(double markov_time) :
        markov_time(markov_time) {
    }

    template<typename G, typename P, typename W>
    double operator ()(G &graph, P &partition, W &weights, double time) {
        cliques::LinearisedInternals internals(graph, weights, partition);
        return (*this)(internals);
    }

    template<typename I>
    double operator ()(I &internals) {
        double q = 1.0 - markov_time;
        int size = internals.comm_w_tot.size();
        cliques::output("time", markov_time, "size", size);


        for (int i = 0; i < size; i++) {
            if (internals.comm_w_tot[i] > 0) {
                q += markov_time * double(
                        internals.comm_w_in[i] / internals.two_m) - ((double(
                        internals.comm_w_tot[i]) / internals.two_m) * (double(
                        internals.comm_w_tot[i]) / internals.two_m));
            }
            cliques::output("q",q);
        }

        return q;
    }

    template<typename I>
    double operator ()(I &internals, int comm_id) {
        double q = -1;
        int i = comm_id;
        if (internals.comm_w_tot[i] > 0) {
//            cliques::output("internals", internals.comm_w_tot[i], internals.two_m, markov_time);
            q = internals.comm_w_tot[i] / internals.two_m * (1.0 - markov_time);
            q += markov_time * double(internals.comm_w_in[i] / internals.two_m)
                    - ((double(internals.comm_w_tot[i]) / internals.two_m)
                        * (double(internals.comm_w_tot[i])/ internals.two_m));
        } else {
            std::cout << "This community does not exist!!!" << std::endl;
        }
//        cliques::output("Q", q);
        return q / internals.comm_w_tot[i] * internals.two_m;
    }

};

/**
 @brief  Functor for finding stability gain for weighted graph
 */
struct linearised_stability_gain_louvain {
    double markov_time;

    linearised_stability_gain_louvain(double mtime) :
        markov_time(mtime) {
    }

    double operator ()(double tot_w_comm, double w_node_to_comm, double two_m,
            double w_deg_node) {
        return (markov_time * w_node_to_comm - tot_w_comm * w_deg_node / two_m);
    }

    template<typename I>
    double operator ()(I &internals, int comm_id_neighbour, int node_id) {
        double tot_w_comm = internals.comm_w_tot[comm_id_neighbour];
        double w_node_to_comm =
                internals.node_weight_to_communities[comm_id_neighbour];
        double w_deg_node = internals.node_to_w[node_id];
        return (markov_time * w_node_to_comm - tot_w_comm * w_deg_node
                / internals.two_m) * 2 / internals.two_m;
    }
};

}
