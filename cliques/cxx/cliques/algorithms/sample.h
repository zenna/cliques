#pragma once

#include <random>

namespace cliques {

template<typename G, typename P>
void uniform_sample(G &graph, int num_nodes, int num_samples, int num_steps_per_sample = 10) {
    typedef typename boost::unordered_set<P, cliques::partition_hash,
            cliques::partition_equal> partition_set;

    int num_sampled = 0;
    int num_steps = 0;
    P current_partition;
    current_partition.initaliase_as_singletons();

    while (true) {
        partition_set neigh_partitions;
        cliques::find_neighbours(graph, current_partition, neigh_partitions);
        int num_current_neighs = neigh_partitions.size();

        while (true) {
            std::uniform_int_distribution randor(0, num_current_neigs);
            P &proposed_partition = randor.find();
            partition_set proposed_neighs;
            cliques::find_neighbours(graph, proposed_partition, proposed_neighs);
            int num_proposed_neighs = proposed_neighs.size();
            float alpha = float(num_current_neighs) / float(num_proposed_neighs);
            std::uniform_real_distribution<double> rand_float(0, 1);
            float rand_num = rand_float();

            if (rand_num > alpha) {
                num_steps++;
                current_partition = proposed_partition;
                break;
            }
        }
        if (num_Steps % num_hops_per_sample == 0) {
            sampled_partitons.insert(current_partition);
            num_sampled++;
        }
        if (num_sampled == num_samples) {
            break;
        }
    }
}

}
