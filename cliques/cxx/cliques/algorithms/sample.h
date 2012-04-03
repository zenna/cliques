#pragma once

#include <random>
#include <math.h>

#include <cliques/algorithms/space.h>

namespace cliques {

/**
 @brief  Use Metropolis Hastings to approximate distribution of stability values
 // TODO make type partition type independent
 */
template<typename G, typename Q, typename M, typename Logger>
void sample_metropolis(G &graph, Q &quality_function,
		double markov_time, int num_samples, int num_steps_per_sample,
		M &sample_to_count, Logger &logger) {
	typedef typename boost::unordered_set<cliques::VectorPartition,
			cliques::partition_hash, cliques::partition_equal> partition_set;

	int num_sampled = 0;
	int num_steps = 0;
	int num_nodes = lemon::countNodes(graph);
	cliques::VectorPartition current_partition(num_nodes);
	current_partition.initialise_as_singletons();

	std::uniform_real_distribution<> real_distribution(0, 1);
	std::mt19937 m_engine;
	std::mt19937 engine; // Mersenne twister MT19937

	while (true) {
		partition_set neigh_partitions;
		cliques::find_neighbours(graph, current_partition, neigh_partitions);
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

			int rand_neigh = distribution(engine);

			// (slow) way to get nth neigh partition
			auto set_itr = neigh_partitions.begin();
			for (int i = 0; i < rand_neigh; ++i) {
				++set_itr;
			}
			cliques::VectorPartition proposed_partition = *set_itr;
			logger.log(proposed_partition);
			partition_set proposed_neighs;

			// Need num neighbours to compute Q(x_t ; x_p)
			cliques::find_neighbours(graph, proposed_partition, proposed_neighs);
			int num_proposed_neighs = proposed_neighs.size();

			//Metropolis-Hastings acceptance alpha
			double prop_energy = quality_function(proposed_partition,
					markov_time);

			double a1 = std::exp(prop_energy - curr_energy);
			double a2 = double(num_current_neighs) / double(
					num_proposed_neighs);

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

// Problems:
////find_neighbours is slow - memoise
//
//template<typename G, typename P, typename Q>
//void stability_dist_sample(G &graph, Q quality_function, double markov_time,
//		int num_samples, int num_steps_per_sample = 1) {
//	typedef typename boost::unordered_set<P, cliques::partition_hash,
//			cliques::partition_equal> partition_set;
//
//	// initialise things and start from singletons
//	int num_sampled = 0;
//	int num_steps = 0;
//	P current_partition;
//	current_partition.initaliase_as_singletons();
//
//	while (true) {
//		// get neighbouring partitions
//		partition_set neigh_partitions;
//		cliques::find_neighbours(graph, current_partition, neigh_partitions);
//		// get energy and neighbours of current state
//		double curr_energy = quality_function(current_partition, markov_time);
//		int num_current_neighs = neigh_partitions.size();
//
//		// Repeat until new proposal is accepted
//		while (true) {
//			// random number generator
//			std::uniform_int_distribution randor(0, num_current_neigs);
//			P &proposed_partition = randor.find();
//
//			double prop_energy = quality_function(proposed_partition,
//					markov_time);
//
//			double a1 = std::exp(prop_energy - curr_energy);
//
//			int num_proposed_neighs = proposed_neighs.size();
////			double a_2 = double(num_proposed_neighs)
////					/ double(num_current_neighs);
//
//			double a_2 = double(num_current_neighs)
//							/ double(num_proposed_neighs);
//
//			double alpha = a1 * a2;
//
//			std::uniform_real_distribution<double> rand_float(0, 1);
//			double rand_num = rand_float();
//
//			if (rand_num < alpha) {
//				num_steps++;
//				current_partition = proposed_partition;
//				break;
//			}
//		}
//		if (num_Steps % num_hops_per_sample == 0) {
//			sampled_partitons.insert(current_partition);
//			num_sampled++;
//		}
//		if (num_sampled == num_samples) {
//			break;
//		}
//	}
//}
//
//


template<typename G, typename P>
void uniform_sample(G &graph, int num_nodes, int num_samples,
		int num_steps_per_sample = 10) {
	typedef typename boost::unordered_set<P, cliques::partition_hash,
			cliques::partition_equal> partition_set;

	// initialise things and start from singletons
	int num_sampled = 0;
	int num_steps = 0;
	P current_partition;
	current_partition.initaliase_as_singletons();

	while (true) {
		// get neighbouring partitions
		partition_set neigh_partitions;
		cliques::find_neighbours(graph, current_partition, neigh_partitions);
		// get number of neigbours
		int num_current_neighs = neigh_partitions.size();

		while (true) {
			// random number generator
			std::uniform_int_distribution randor(0, num_current_neigs);
			P &proposed_partition = randor.find();
			partition_set proposed_neighs;
			cliques::find_neighbours(graph, proposed_partition, proposed_neighs);
			int num_proposed_neighs = proposed_neighs.size();
			float alpha = float(num_current_neighs)
					/ float(num_proposed_neighs);
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
