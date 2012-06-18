    /* Copyright (c) Z Tavares, M Schaub - zennatavares@gmail.com, 2010-2012
 Stochastic hill climbing optimisation functions */
#pragma once
#include <functional>

#include <cliques/helpers/stochastic.h>

namespace cliques {
/**
 @brief  Strongly typed enum for direction of optimisation (min/max)
 */
enum class Direction {
    ASCENT, DESCENT
};

/**
 @brief  Given a configuration, climb stochastically and monotonically until local optimum is reached
 @tparam C                          Configuration type (e.g. partition, program)
 @tparam cC                         Container for functors
 @tparam QF                         Quality functor type
 @param[in]  initial_configuration  Configuration - i.e. solution  (e.g. partition, program)
 @param[in]  find_neighbours        function giving local neighbourhood of configuration
 @param[in]  direction              Ascent (maximisation) or Descent (minimisation)
 @param[in]  compute_quality        Computes the quality (or cost) of a given configuration
 @param[out] local_optimum          Local optimum found through schocastic hill climbing
 */
template<typename C, typename cC, typename RNG>
C stochastic_monotonic_climb(
        C const initial_configuration,
        std::function<cC (C)> const find_neighbours,
        cliques::Direction const direction,
        std::function<double (C)> const compute_quality,
        RNG &prng_engine) {

    bool am_at_local_optimum = false; // stop when this is true
    int max_num_steps = 10000;
    C current_configuration = initial_configuration;

    // Make stochastic local move until local optima or max_num_steps reached
    for (int i = 0; i < max_num_steps && am_at_local_optimum == false; ++i) {

        // First find differences in quality between current config and all neighbours
        double current_config_quality = compute_quality(current_configuration);
        std::vector<double> quality_diffs;

        // cliques::output("my quality:",current_config_quality);
        // cliques::print_partition_line(current_configuration);

        cC neighbours = find_neighbours(current_configuration);
        
        am_at_local_optimum = true;
        for (C const &neighbour : neighbours) {
            // cliques::output("neighbour quality:",compute_quality(neighbour));
            double quality_diff;
            if (direction == cliques::Direction::ASCENT) {
                quality_diff =  compute_quality(neighbour) - current_config_quality;
            }
            else if (direction == cliques::Direction::DESCENT) {
                quality_diff = current_config_quality - compute_quality(neighbour) ;
            }
            // cliques::output("quality_diff", quality_diff);

            am_at_local_optimum = am_at_local_optimum && (quality_diff <= 0);
            
            // We don't want to sample worse neighbours, so set them to zero
            quality_diff = quality_diff > 0 ? quality_diff : 0;
            quality_diffs.push_back(quality_diff);
        }

        if (am_at_local_optimum == true) {
            cliques::output("got to local max in steps:", i+1);
        }
        else if (i == max_num_steps-1) {
            cliques::output("didn't reach local maxima");
        }

        // Then sample a neighbour from multinomial distribution of qualities
        // Return if we've found a local optimum
        if (am_at_local_optimum == false) {
            // cliques::output("moving to:");
            int chosen_index = weighted_sample(quality_diffs, prng_engine);
            typename cC::iterator it = neighbours.begin();
            std::advance(it, chosen_index);
            // cliques::output("diff is:", compute_quality(*it) - current_config_quality);
            current_configuration = *it;
            // cliques::print_partition_line(current_configuration);
        }

        cliques::output("level: ",i, " quality:",compute_quality(current_configuration));
    }

    return current_configuration;
}

}
