#ifndef CLIQUES_MODULE_H
#define CLIQUES_MODULE_H

#include <limits>
#include <vector>
#include <lemon/concepts/graph.h>

#include <cliques/structures/partition.h>
#include <cliques/structures/common.h>

namespace cliques {

namespace _luv {

inline void insert_node_bookkeeping(int node, int comm, double)

// TODO insert node, remove node functions...


}
/**
 @brief  Louvain method - greedy algorithm to find community structure of a network.

 This is a fast algorithm to find communities within a network
 ref: Fast Unfolding Of Communities in large networks, Blondel et al. 2008

 @param[in]  my_graph     graph to find partition of
 @param[in]  quality_function     partition quality function object
 */
template<typename P, typename T, typename QF>
P find_optimal_partition_louvain(T &graph, QF quality_function) {
	// Create singleton partition from graph
	P partition;
	for (typename T::NodeIt n(graph); n != lemon::INVALID; ++n) {
		partition.add_element(graph.id(n));
	}

	for (typename T::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		int best_set = 0;
		float best_stability = -std::numeric_limits<float>::max();

		for (typename T::IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
			typename T::Node n2 = graph.oppositeNode(n1, e);
			if (n1 != n2) {
				//std::cout << "trying " << graph.id(n1) << " " << graph.id(n2) << std::endl;

				partition.union_sets(graph.id(n1), graph.id(n2));

				std::vector<float> new_stability;
				quality_function(graph, partition, new_stability);
				partition.undo_last_union();

				//std::cout << "new stab "<< new_stability[0] << "best stab: " << best_stability << std::endl;
				if (new_stability[0] > best_stability) {
					best_stability = new_stability[0];
					best_set = graph.id(n2);
				}
			}
		}
		std::vector<float> old_stability;
		quality_function(graph, partition, old_stability);
		//std::cout << "original stab "<< old_stability[0] << std::endl;
		if (best_stability > old_stability[0]) {
			std::cout << "joining " << graph.id(n1) << " " << best_set
					<< std::endl;
			partition.union_sets(graph.id(n1), best_set);
		}
	}

	return partition;
}

/* Implementation for Luvain algorithm for "modularity-like" quality functions, where the change
 * in quality can be evaluated efficiently.
 */

template<typename P, typename T, typename W, typename QF, typename QFDIFF>
P find_optimal_partition_louvain(T &graph, W &weights, QF quality_function,
		QFDIFF quality_function_diff) {

	// Start: Create singleton partition from graph
	P partition;
	for (typename T::NodeIt n(graph); n != lemon::INVALID; ++n) {
		partition.add_element(graph.id(n));
	}

	// some bookkeeping	to keep matters simple
	// total weight of the links
	two_m = find_total_weight(&graph, &weights);
	// number_of_nodes
	number_of nodes = lemon::countNodes(graph);
	// mapping form node to node weight
	node_to_w = lemon::rangeMap(number_of_nodes);
	// mapping from node number to community number
	node_to_comm = lemon::RangeMap<int>(number_of_nodes, 0);
	// mapping #community to total weight of community
	comm_w_tot = lemon::RangeMap<int>(number_of_nodes, 0);
	// mapping #community to weight within community
	comm_w_in = lemon::RangeMap<int>(number_of_nodes, 0);

	for (int i = 0; i < number_of_nodes; i++) {
		temp_node = graph.nodeFromId(i);
		comm_w_tot[i] = node_to_w[i] = find_weighted_degree(graph, weights,
				temp_node);
		node_to_comm[i] = i;
		comm_w_in[i] = find_
	}

	for (typename T::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		int best_set = 0;
		float best_stability = -std::numeric_limits<float>::max();

		for (typename T::IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
			typename T::Node n2 = graph.oppositeNode(n1, e);
			if (n1 != n2) {
				//std::cout << "trying " << graph.id(n1) << " " << graph.id(n2) << std::endl;
				partition.union_sets(graph.id(n1), graph.id(n2));

				std::vector<float> new_stability;
				quality_function(graph, partition, new_stability);
				partition.undo_last_union();

				//std::cout << "new stab "<< new_stability[0] << "best stab: " << best_stability << std::endl;
				if (new_stability[0] > best_stability) {
					best_stability = new_stability[0];
					best_set = graph.id(n2);
				}
			}
		}
		std::vector<float> old_stability;
		quality_function(graph, partition, old_stability);
		//std::cout << "original stab "<< old_stability[0] << std::endl;
		if (best_stability > old_stability[0]) {
			std::cout << "joining " << graph.id(n1) << " " << best_set
					<< std::endl;
			partition.union_sets(graph.id(n1), best_set);
		}
	}

	return partition;
}

}

#endif
