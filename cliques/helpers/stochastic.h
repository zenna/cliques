/* Copyright (c) Z Tavares, M Schaub - zennatavares@gmail.com, 2010-2012
 Helper functions for sampling and stochastic operations */
#pragma once
#include <cliques/helpers/lists.h>

namespace cliques {

/**
 @brief  Samples from multinomial distribution
 @tparam C                          Configuration type (e.g. partition, program)
 @tparam R                          PRNG engine
 @param[in]  weighted_list          list of weights, e.g. [9.3,0.1,100]
 @param[in]  prng_engine            PRNG engine - stores state of RNG, e.g. std::mt19937
 @param[out] index_of_sample        Vector containing the sampled partitions
 */
template<typename N, typename RNG>
int weighted_sample(std::vector<N> weighted_list, RNG &prng_engine) {
    std::uniform_real_distribution<> real_distribution(0, 1);

    // indexed_array has each element as sum of previous elements in weighted_list
    // i.e. on line from 0 to total weight, assign each individual weight adequate fraction
    double total_weight = sum(weighted_list);
    std::vector<N> indexed_array;
    N total = 0;
    for (N &weight : weighted_list) {
        total += weight;
        indexed_array.push_back(total);
    }

    // Then sample random point on this line, and see corresponding weight
    double rand_real = real_distribution(prng_engine) * total_weight;

    int rand_index = 0;
    for (N &list_element : indexed_array) {
        if (rand_real <= list_element) {
            return rand_index;
        }
    rand_index += 1;
    }
}

}
