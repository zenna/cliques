/* Copyright (c) Zenna Tavares - zennatavares@gmail.com */
#pragma once
#include "cliques/algorithms/louvain.h"
#include "cliques/algorithms/varofinf.h"
#include "cliques/quality_functions/stability.h"
#include "cliques/structures/vector_partition.h"

namespace cliques {

template<typename G, typename M>
double variance_of_vi(G &graph, M &weights, double markov_time,
        int num_iterations) {
    typedef cliques::VectorPartition VectorPartition;
    std::vector<VectorPartition> all_optima;
    VectorPartition singleton_partition(lemon::countNodes(graph));
    singleton_partition.initialise_as_singletons();

    // Find optima with louvain
    for (int i = 0; i < num_iterations; ++i) {
        double current_markov_time = markov_time;
        std::vector<double> markov_times(1, markov_time);
        cliques::NoLogging no_logging;
        std::vector<VectorPartition> optimal_partitions;

        cliques::find_optimal_partition_louvain_with_gain(
                graph,
                weights,
                cliques::find_weighted_linearised_stability(markov_times),
                cliques::linearised_stability_gain_louvain(current_markov_time),
                singleton_partition, optimal_partitions, no_logging);

        VectorPartition optima = optimal_partitions.back();
        all_optima.push_back(optima);
    }

    std::vector<double> var_of_infs;

    // Compute all pairwise var_of_infs
    double total = 0.0;
    int num_elements = (all_optima.size() * (all_optima.size() - 1)) / 2;
    for (auto p1 = all_optima.begin(); p1 != all_optima.end(); ++p1) {
        for (auto p2 = p1; ++p2 != all_optima.end();) {
            double vi = find_variation_of_information(*p1, *p2);
            var_of_infs.push_back(vi);
            total += vi / num_elements;
        }
    }

    return total;
}

template<typename G, typename M>
double basin_size(G &graph, M &weights, double markov_time) {

    return 0.3;
}

}
