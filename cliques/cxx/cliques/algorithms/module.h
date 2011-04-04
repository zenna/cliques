#ifndef CLIQUES_MODULE_H
#define CLIQUES_MODULE_H

#include <limits>
#include <vector>
#include <lemon/concepts/graph.h>

#include <cliques/structures/partition.h>
#include <cliques/structures/common.h>

namespace cliques
{
/**
@brief  Louvain method - greedy algorithm to find community structure of a network.

This is a fast algorithm to find communities within a network
ref: Fast Unfolding Of Communities in large networks, Blondel et al. 2008

@param[in]  my_graph     graph to find partition of
@param[in]  quality_function     partition quality function object
*/
template <typename P, typename T, typename QF>
P find_optimal_partition_louvain(T &graph,
								QF quality_function)
{
	// Create singleton partition from graph
	P partition;
	for (typename T::NodeIt n(graph); n!= lemon::INVALID; ++n) {
		partition.add_element(graph.id(n));
	}

	for (typename T::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		int best_set = 0;
		float best_stability = -std::numeric_limits<float>::max();

		for (typename T::IncEdgeIt e(graph,n1); e != lemon::INVALID; ++e) {
			typename T::Node n2 = graph.oppositeNode (n1, e);
			if (n1 != n2) {
				//std::cout << "trying " << graph.id(n1) << " " << graph.id(n2) << std::endl;

				partition.union_sets(graph.id(n1),graph.id(n2));

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
			std::cout << "joining " << graph.id(n1) << " " << best_set << std::endl;
			partition.union_sets(graph.id(n1), best_set);
		}
	}

	return partition;
}
}

#endif
