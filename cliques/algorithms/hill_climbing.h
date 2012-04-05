/* Copyright (c) Z Tavares, M Schaub - zennatavares@gmail.com, 2010-2012
 Stochastic hill climbing optimisation functions */
#pragma once
#include <functional>

#include <cliques/helpers/stochastic.h>

namespace cliques {
/**
 @brief  Strongly typed enum for direction of stochastic ascent
 */
enum class Direction {
    ASCENT, DESCENT
};

/**
 @brief  Given a configuration, climb stochastically and monotonically until local optimum is reached
 @tparam C                          Configuration type (e.g. partition, program)
 @tparam QF                         Quality functor type
 @param[in]  initial_configuration  Configuration - i.e. solution  (e.g. partition, program)
 @param[in]  find_neighbours        function giving local neighbourhood of configuration
 @param[in]  direction              Ascent (maximisation) or Descent (minimisation)
 @param[in]  compute_quality        Computes the quality (or cost) of a given configuration
 @param[out] local_optimum           Vector containing the sampled partitions
 */
template<typename C>
C stochastic_monotonic_climb(
        C initial_configuration,
        std::function<std::vector<C> (C)> find_neighbours,
        cliques::Direction direction,
        std::function<double (C)> compute_quality) {

    int max_num_steps = 10000;
    C current_configuration = initial_configuration;
    std::mt19937 m_engine;

    // Make stochastic local move until local optima or max_num_steps reached
    for (int i = 0; i < max_num_steps; ++i) {

        // First find differences in quality between current config and all neighbours
        double current_config_quality = compute_quality(initial_configuration);
        std::vector<double> quality_diffs;
        std::vector<C> neighbours = find_neighbours(current_configuration);
        for (C &neighbour : neighbours) {
            double quality_diff;
            if (direction == cliques::Direction::ASCENT) {
                quality_diff =  compute_quality(current_configuration) - current_config_quality;
            }
            else if (direction == cliques::Direction::DESCENT) {
                quality_diff = current_config_quality - compute_quality(current_configuration);
            }

            // By setting weight of negative quality diffs to 0, we enforce monotonic moves
            // Since the weighted sampler will not chose these
            quality_diff = quality_diff > 0 ? quality_diff : 0;
            quality_diffs.push_back(quality_diff);
        }

        // Then sample a neighbour from multinomial distribution of qualities
        // Return if we've found a local optimum
        if (quality_diffs.size() > 0) {
            int chosen_index = weighted_sample(quality_diffs, m_engine);
            auto current_configuration = neighbours[chosen_index];
        } else {
            return current_configuration;
        }
    }
}

}
