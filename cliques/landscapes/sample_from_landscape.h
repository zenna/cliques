#pragma once

#include <random>
#include <math.h>

#include <cliques/landscapes/landscape_space.h>
#include <lemon/concepts/graph.h>

namespace clq {

/**
 @brief  Use Metropolis Hastings to approximate distribution of stability values
 // TODO make type partition type independent
 */
template<typename G, typename Q, typename M, typename Logger>
void sample_metropolis(G &graph, Q &quality_function, double markov_time,
        int num_samples, int num_steps_per_sample, M &sample_to_count,
        Logger &logger) {
    typedef typename std::unordered_set<clq::VectorPartition,
            clq::partition_hash, clq::partition_equal> partition_set;

    int num_sampled = 0;
    int num_steps = 0;
    int num_nodes = lemon::countNodes(graph);
    clq::VectorPartition current_partition(num_nodes);
    current_partition.initialise_as_singletons();

    std::uniform_real_distribution<> real_distribution(0, 1);
    std::mt19937 m_engine; // Mersenne twister MT19937

    while (true) {
        partition_set neigh_partitions = clq::find_neighbours(graph, current_partition);
        int num_current_neighs = neigh_partitions.size();
        std::uniform_int_distribution<int> distribution(0,
                num_current_neighs - 1);

        double curr_energy = quality_function(current_partition, markov_time);

        while (true) {
            // Record sample if at multiple of num_steps_per_sample
            if (num_steps % num_steps_per_sample == 0) {
                current_partition.normalise_ids();
                sample_to_count[current_partition]++;
                logger.log(current_partition);
                num_sampled++;
            }

            int rand_neigh = distribution(m_engine);

            // (slow) way to get nth neigh partition
            auto set_itr = neigh_partitions.begin();
            for (int i = 0; i < rand_neigh; ++i) {
                ++set_itr;
            }
            clq::VectorPartition proposed_partition = *set_itr;
            logger.log(proposed_partition);
            partition_set proposed_neighs;

            // Need num neighbours to compute Q(x_t ; x_p)
            clq::find_neighbours(graph, proposed_partition);
            int num_proposed_neighs = proposed_neighs.size();

            //Metropolis-Hastings acceptance alpha
            double prop_energy = quality_function(proposed_partition,
                    markov_time);

            double a1 = std::exp(prop_energy - curr_energy);
            double a2 = double(num_current_neighs)
                    / double(num_proposed_neighs);

            double alpha = a1 * a2;

            double rand_uniform_num = real_distribution(m_engine);

            // Increment time and do move to new state if accepted
            num_steps++;
            if (rand_uniform_num < alpha) {
                current_partition = proposed_partition;
                break;
            }
        }

        if (num_sampled >= num_samples) {
            break;
        }
    }
}

/**
 @brief  Sample uniformly from the set of all allowed partitions using Metropolis Hastings MCMC method.
 @tparam P                          Partition type
 @tparam G                          Graph type
 @param[in]  graph     				Graph to partition (not the landscape graph)
 @param[in]  num_samples 			Number of samples that should be drawn
 @param[in]  num_steps_per_sample	Number of steps before sample is taken as output
 @param[out] output_partitions		Vector containing the sampled partitions
 */
template<typename P, typename G>
std::vector<P> uniform_sample(G &graph, int num_samples,
        int num_steps_per_sample = 10) {

    // define partition_set as unordered_set of partitions
    typedef typename std::unordered_set<P, clq::partition_hash,
            clq::partition_equal> partition_set;

    std::vector<P> sampled_partitions;

    // initialise variables and start from singletons
    int num_sampled = 0;
    int num_steps_taken = 0;
    P current_partition(lemon::countNodes(graph));
    current_partition.initialise_as_singletons();

    // initialise random number generator
    std::uniform_real_distribution<> real_distribution(0, 1);
    std::mt19937 m_engine; // Mersenne twister MT19937

    // initialise and find neighbouring partitions
    partition_set neigh_partitions = clq::find_neighbours(graph, current_partition);

    // get number of neighbours of current partition
    int num_current_neighs = neigh_partitions.size();

    // start "random walk" on the partition landscape until enough samples are collected
    while (num_sampled < num_samples) {
        // initialise uniform distribution over neighbours and draw random neighbour
        std::uniform_int_distribution<int> distribution(0,
                num_current_neighs - 1);
        int rand_neigh = distribution(m_engine);

        // get the randomly picked partition (slow)
        auto set_itr = neigh_partitions.begin();
        for (int i = 0; i < rand_neigh; ++i) {
            ++set_itr;
        }
        P proposed_partition = *set_itr;

        // find the neighbours of the proposed partition, to compute acceptance ratio etc.
        partition_set proposed_neighs = clq::find_neighbours(graph, proposed_partition);
        int num_proposed_neighs = proposed_neighs.size(); // numbber of proposed neighbours

        //Metropolis-Hastings acceptance probability
        double alpha = double(num_current_neighs) / double(num_proposed_neighs);

        // Move to new state if accepted and update statistics of the state
        if (alpha >= 1 || real_distribution(m_engine) < alpha) {
            current_partition = proposed_partition;
            num_current_neighs = num_proposed_neighs;
            neigh_partitions.swap(proposed_neighs); //TODO: check if there is not a faste way of doing this, direct initialisation did not work however
        }

        // Record sample if at multiple of num_steps_per_sample
        if (num_steps_taken % num_steps_per_sample == 0) {
            current_partition.normalise_ids(); // TODO check if this is needed
            sampled_partitions.push_back(current_partition);
            num_sampled++;
        }
        num_steps_taken++;
    }
    return sampled_partitions;
}

}
