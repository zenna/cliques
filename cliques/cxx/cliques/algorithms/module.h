/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_MODULE_H
#define CLIQUES_MODULE_H

#include <limits>
#include <vector>
#include <lemon/concepts/graph.h>

#include <cliques/structures/partition.h>
#include <cliques/structures/common.h>

namespace cliques {

namespace _luv {


inline void insert_node_bookkeeping(int node, int comm, double){

}

// TODO insert node, remove node functions...


}
/**
@brief  Louvain method - greedy algorithm to find community structure of a network.

This is a fast algorithm to find communities within a network
ref: Fast Unfolding Of Communities in large networks, Blondel et al. 2008

@param[in]  graph     graph to find partition of
@param[in]  compute_quality     partition quality function object
@param[out]  std::vector<partitions>     optimal partitions, last in vector is overall best
*/
template <typename P, typename M, typename G, typename QF>
P find_optimal_partition_louvain(G &graph,
								QF compute_quality,
								M &weights,
								P &optimal_partitions)
{
	bool did_quality_increase = false;
	// Create singleton partition from graph
	P partition;
	for (typename G::NodeIt n(graph); n!= lemon::INVALID; ++n) {
		partition.add_element(graph.id(n));
	}

	for (typename G::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		int best_set = 0;
		float best_stability = -std::numeric_limits<float>::max();

		for (typename G::IncEdgeIt e(graph,n1); e != lemon::INVALID; ++e) {
			typename G::Node n2 = graph.oppositeNode (n1, e);
			if (n1 != n2) {
				//std::cout << "trying " << graph.id(n1) << " " << graph.id(n2) << std::endl;

				partition.union_sets(graph.id(n1), graph.id(n2));

				std::vector<float> new_stability;
				compute_quality(graph, partition, weights, new_stability);
				partition.undo_last_union();

				//std::cout << "new stab "<< new_stability[0] << "best stab: " << best_stability << std::endl;
				if (new_stability[0] > best_stability) {
					best_stability = new_stability[0];
					best_set = graph.id(n2);
				}
			}
		}
		std::vector<float> old_stability;
		compute_quality(graph, partition, old_stability);
		//std::cout << "original stab "<< old_stability[0] << std::endl;
		if (best_stability > old_stability[0]) {
			std::cout << "joining " << graph.id(n1) << " " << best_set
					<< std::endl;
			partition.union_sets(graph.id(n1), best_set);
			did_quality_increase = true;
		}
	}

	// The second phase of the algorithm consists in building a new network
	// whose nodes are now the communities found during the first phase.
	// To do so, the weights of the links between the new nodes are given by
	// the sum of the weight of the links between nodes in the corresponding
	// two communities.
	// Links between nodes of the same community lead to self-loops for this
	// community in the new network. Once this second phase is completed,
	// it is then possible to reapply the first phase of the algorithm
	// to the resulting weighted network and to iterate.
	if (did_quality_increase == true) {
		M reduced_weight_map;
		// Compile P and graph back into normal partition
		// Pushback P into list of partitions
		G reduced_graph;

		for (typename G::NodeIt node(graph); node != lemon::INVALID; ++node) {
			reduced_graph.addNode();
		}
		//Need a map set_id > node_in_reduced_graph_

		// Find between community total weights by checking
		// Edges within graph
		for (typename G::EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
			int comm_of_node_u = partition.find_set(graph.id(graph.u(edge)));
			int comm_of_node_v = partition.find_set(graph.id(graph.v(edge)));

			if (comm_of_node_u != comm_of_node_v) {
				float weight = weights[edge];
				typename G::Edge edge_in_reduced_graph = lemon::findEdge(reduced_graph, graph.nodeFromId(comm_of_node_u),
						graph.nodeFromId(comm_of_node_u));
				reduced_weight_map[edge_in_reduced_graph] += weight;
			}
		}

		return cliques::find_optimal_partition_louvain<P>(graph, reduced_weight_map,
				compute_quality, optimal_partitions);
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
		comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node);
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
