/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010- 2011 */
#pragma once

#include <limits>
#include <vector>
#include <set>

#include <lemon/concepts/graph.h>

#include <cliques/algorithms/internals/internals.h>

namespace cliques {
/**
 @brief  Finds the maxima of a landscape graph through steepest ascent

 This is an exhaustive search to find maximal points on a landscape
 when a landscape is represented as a graph

 @param[in] graph the input graph
 @param[in] stabilities pointer (array) of node 'energies'
 @tparam G graph type
 */
template<typename G>
std::set<int> find_maxima(G &graph, double *stabilities) {
    std::set<int> maxima;
    int num_iterations = 0;
    for (typename G::NodeIt n(graph); n != lemon::INVALID; ++n) {
        // Check part_itr hasn't already been visited
        if (num_iterations % 100 == 0) {
            //std::cout << "the number of iterations is" << num_iterations << std::endl;
            //std::cout << "so far number of maxima is " << maxima.size() << std::endl;
            //std::cout << "size of visited is " << visited.size() << std::endl;
            //std::cout << "num skipped " << num_skipped << std::endl;
        }
        typename G::Node best_neighbour = n;
        while (1) {
            bool has_improved = false;
            double best_score = stabilities[graph.id(best_neighbour)];
            for (typename G::OutArcIt a(graph, best_neighbour); a
                    != lemon::INVALID; ++a) {
                if (best_score < stabilities[graph.id(graph.target(a))]) {
                    best_score = stabilities[graph.id(graph.target(a))];
                    best_neighbour = graph.target(a);
                    has_improved = true;
                }
            }

            if (has_improved == false) {
                maxima.insert(graph.id(best_neighbour));
                break;
            }
        }
        num_iterations++;
    }
    return maxima;
}

/**
 @brief  Samples the maxima of a landscape graph through steepest ascent

 This is an sampled search to find maximal points on a landscape
 */
template<typename T, typename W, typename QF, typename QFDIFF, typename P, typename Logger>
void sample_maxima(T &graph, W &weights, QF compute_quality,
        QFDIFF compute_quality_diff, boost::unordered_set<P,
                cliques::partition_hash, cliques::partition_equal> &maxima,
        boost::unordered_set<P, cliques::partition_hash,
                cliques::partition_equal> &sampled_partitions, Logger log) {

    typedef typename boost::unordered_set<P, cliques::partition_hash,
            cliques::partition_equal> partition_set;

    auto internals = cliques::gen(compute_quality, graph, weights);
    P partition(lemon::countNodes(graph));
    partition.initialise_as_singletons();

    for (auto set_itr = sampled_partitions.begin(); set_itr
            != sampled_partitions.end(); ++set_itr) {
        auto best_neighbour = *set_itr;
        while (true) {
            bool has_improved = false;
            partition_set neighs;
            cliques::find_neighbours(graph, best_neighbour, neighs);
            auto internals = cliques::gen(compute_quality, graph, weights, best_neighbour);
            double best_quality = compute_quality(internals);
            for (auto neigh_itr = neighs.begin(); neigh_itr != neighs.end(); ++neigh_itr) {
                auto neigh_internals = cliques::gen(compute_quality, graph, weights, *neigh_itr);
                //Internals neigh_internals(graph, weights, *neigh_itr);
                double neigh_quality = compute_quality(neigh_internals);
                if (best_quality < neigh_quality) {
                    best_quality = neigh_quality;
                    best_neighbour = *neigh_itr;
                    has_improved = true;
                }
            }
            if (has_improved == false) {
            	best_neighbour.normalise_ids();
                maxima.insert(best_neighbour);
                break;
            }
        }
    }
}

}
