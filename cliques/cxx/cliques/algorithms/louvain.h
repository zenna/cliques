/* Copyright (c) Z Tavares, M Schaub - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_LOUVAIN_H
#define CLIQUES_LOUVAIN_H

#include <limits>
#include <vector>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/graphhelpers.h>
#include <cliques/structures/partition.h>
#include <cliques/structures/common.h>
#include <iostream>
#include <ctime>
#include <boost/random.hpp>

namespace cliques {

/**
 @brief  Louvain method - greedy algorithm to find community structure of a network.

 This is a fast algorithm to find communities within a network
 ref: Fast Unfolding Of Communities in large networks, Blondel et al. 2008

 The second phase of the algorithm consists in building a new network
 whose nodes are now the communities found during the first phase.
 To do so, the weights of the links between the new nodes are given by
 the sum of the weight of the links between nodes in the corresponding
 two communities.
 Links between nodes of the same community lead to self-loops for this
 community in the new network. Once this second phase is completed,
 it is then possible to reapply the first phase of the algorithm
 to the resulting weighted network and to iterate.

 @param[in]  graph     graph to find partition of
 @param[in]  compute_quality     partition quality function object
 @param[out]  std::vector<partitions>     optimal partitions, last in vector is overall best
 */
template<typename P, typename M, typename G, typename QF>
P find_optimal_partition_louvain(G &graph, M &weights, QF compute_quality,
		std::vector<P> &optimal_partitions) {

	bool did_quality_increase = false;
	// Create singleton partition from graph
	P partition;
	for (typename G::NodeIt n(graph); n != lemon::INVALID; ++n) {
		partition.add_element(graph.id(n));
	}

	for (typename G::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		int best_set = 0;
		float best_stability = -std::numeric_limits<float>::max();

		for (typename G::IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
			typename G::Node n2 = graph.oppositeNode(n1, e);
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
		compute_quality(graph, partition, weights, old_stability);
		//std::cout << "original stab "<< old_stability[0] << std::endl;
		if (best_stability > old_stability[0]) {
			std::cout << "joining " << graph.id(n1) << " " << best_set
					<< std::endl;
			partition.union_sets(graph.id(n1), best_set);
			did_quality_increase = true;
		}
	}

	if (did_quality_increase == true) {
		M reduced_weight_map(graph);
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
				typename G::Edge edge_in_reduced_graph = lemon::findEdge(
						reduced_graph, graph.nodeFromId(comm_of_node_u),
						graph.nodeFromId(comm_of_node_u));
				reduced_weight_map[edge_in_reduced_graph] += weight;
			}
		}

		return cliques::find_optimal_partition_louvain<P>(graph,
				reduced_weight_map, compute_quality, optimal_partitions);
	}

	return partition;
}

/**
 @brief  isolate a node into its singleton set & update internals
 */
template<typename G, typename M, typename I, typename P>
void isolate_and_update_internals(G &graph, M &weights, typename G::Node node,
		I &internals, P &partition) {
	int node_id = graph.id(node);
	int comm_id = partition.find_set(node_id);
	// get weights from node to each community
	internals.node_weight_to_communities
			= cliques::find_weight_node_to_communities(graph, partition,
					weights, node);
	internals.comm_w_tot[comm_id] -= internals.node_to_w[node_id];
	internals.comm_w_in[comm_id] -= 2
			* internals.node_weight_to_communities[comm_id]
			+ find_weight_selfloops(graph, weights, node);
	partition.isolate_node(node_id);
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename G, typename M, typename I, typename P>
void insert_and_update_internals(G &graph, M &weights, typename G::Node node,
		I &internals, P &partition, int best_comm) {
	int node_id = graph.id(node);
	// insert node to partition/bookkeeping
	//			std::cout << node_id << comm_id << best_comm << std::endl;
	internals.comm_w_tot[best_comm] += internals.node_to_w[node_id];
	//			std::cout << comm_w_tot[best_comm] << std::endl;
	internals.comm_w_in[best_comm] += 2
			* internals.node_weight_to_communities[best_comm]
			+ find_weight_selfloops(graph, weights, node);
	//			std::cout << comm_w_in[best_comm] << std::endl;
	//			std::cout << "was here" << std::endl;
	partition.add_node_to_set(node_id, best_comm);
}

struct Internals {
	typedef lemon::RangeMap<double> range_map;
	unsigned int num_nodes;
	double two_m;
	range_map node_to_w;
	range_map comm_w_tot;
	range_map comm_w_in;
	std::map<int, double> node_weight_to_communities;

	template<typename G, typename M>
	Internals(G &graph, M &weights) :
		num_nodes(lemon::countNodes(graph)), node_to_w(num_nodes, 0),
				comm_w_tot(num_nodes, 0), comm_w_in(num_nodes, 0) {
		two_m = 2 * find_total_weight(graph, weights);
		for (unsigned int i = 0; i < num_nodes; ++i) {
			typename G::Node temp_node = graph.nodeFromId(i);
			comm_w_tot[i] = node_to_w[i] = find_weighted_degree(graph, weights,
					temp_node);
			comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node);
		}
	}

	template<typename G, typename M, typename P>
	Internals(G &graph, M &weights, P &partition) :
		num_nodes(lemon::countNodes(graph)), node_to_w(num_nodes, 0),
				comm_w_tot(num_nodes, 0), comm_w_in(num_nodes, 0) {
		two_m = 2 * find_total_weight(graph, weights);

		typedef typename G::EdgeIt EdgeIt;
		// find internal statistics based on graph, weights and partitions
		// consider all edges
		for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
			int node_u_id = graph.id(graph.u(edge));
			int node_v_id = graph.id(graph.v(edge));

			// this is to distinguish within community weight with total weight
			int comm_of_node_u = partition.find_set(node_u_id);
			int comm_of_node_v = partition.find_set(node_v_id);

			// weight of edge
			double weight = weights[edge];

			// add weight to node weight
			node_to_w[node_u_id] += weight;
			node_to_w[node_v_id] += weight;
			// add weight to total weight of community
			comm_w_tot[comm_of_node_u] += weight;
			comm_w_tot[comm_of_node_v] += weight;
			if (comm_of_node_u == comm_of_node_v) {
				// in case the weight stems from within the community add to internal weights
				comm_w_in[comm_of_node_u] += 2 * weight;
			}
		}
	}
};
/**
 @brief  Louvain method - greedy algorithm to find community structure of a network.

 This specialisation is for efficiency.

 @param[in]  graph     graph to find partition of
 @param[in]  compute_quality     partition quality function object
 @param[out]  std::vector<partitions>     optimal partitions, last in vector is overall best
 */
/*
 * add static bool has_gain_function = true
 * still need to select gain function
 */
template<typename P, typename T, typename W, typename QF, typename QFDIFF>
double find_optimal_partition_louvain_with_gain(T &graph, W &weights,
		QF compute_quality, QFDIFF compute_quality_diff,
		std::vector<P> &optimal_partitions) {

	typedef typename T::Node Node;
	typedef typename T::Edge Edge;
	typedef typename T::NodeIt NodeIt;
	typedef typename T::EdgeIt EdgeIt;
	typedef typename T::IncEdgeIt IncEdgeIt;
	typedef lemon::RangeMap<double> range_map;

	Internals internals(graph, weights);
	P partition(lemon::countNodes(graph));
	partition.initialise_as_singletons();
	double minimum_improve = 0.000001;
	double current_quality = compute_quality(internals);
	bool one_level_end = false;
	double old_quality = current_quality;
	bool did_nodes_move = false;
	//one level louvain

	// Randomise the looping over nodes
	// initialise random number generator
	srand(std::time(0));
	// create vector to shuffle
	std::vector<Node> nodes_ordered_randomly;
	for (NodeIt temp_node(graph); temp_node != lemon::INVALID; ++temp_node) {
		nodes_ordered_randomly.push_back(temp_node);
	}
	// actual shuffling
	std::random_shuffle(nodes_ordered_randomly.begin(),
			nodes_ordered_randomly.end());

	do {
		did_nodes_move = false;
		old_quality = current_quality;

		// loop over all nodes in random order
		typename std::vector<Node>::iterator n1_it;
		for (n1_it = nodes_ordered_randomly.begin(); n1_it
				!= nodes_ordered_randomly.end(); ++n1_it) {
			// get node id and comm id
			Node n1 = *n1_it;
			unsigned int node_id = graph.id(n1);
			unsigned int comm_id = partition.find_set(node_id);
			isolate_and_update_internals(graph, weights, n1, internals,
					partition);
			//default option for re-inclusion of node
			unsigned int best_comm = comm_id;
			double best_gain = 0;

			// TODO check if better to loop over all neighbouring  communities,
			// find all neighbouring communities
			// loop over all neighbouring nodes
			for (IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
				Node n2 = graph.oppositeNode(n1, e);
				// get neighbour node id and neighbour community id
				unsigned int node_id_neighbour = graph.id(n2);
				if (node_id != node_id_neighbour) {
					unsigned int comm_id_neighbour = partition.find_set(
							node_id_neighbour);
					double gain = compute_quality_diff(internals,
							comm_id_neighbour, node_id);
					if (gain > best_gain) {
						best_comm = comm_id_neighbour;
						best_gain = gain;
						// avoid not necessary movements, place node in old community if possible
					} else if (gain == best_gain && comm_id
							== comm_id_neighbour) {
						best_comm = comm_id;
					}
				}
			}

			insert_and_update_internals(graph, weights, n1, internals,
					partition, best_comm);
			if (best_comm != comm_id) {
				did_nodes_move = true;
			}
		}
		if (did_nodes_move == true) {
			current_quality = compute_quality(internals);
			one_level_end = true;
		}

		// If there was a movement of the nodes AND quality increases make another run
	} while ((current_quality - old_quality) > minimum_improve);

	// Start Second phase - create reduced graph with self loops
	partition.normalise_ids();
	if (one_level_end == true) {
		// Compile P into original partition size
		short hierarchy = optimal_partitions.size();
		if (hierarchy == 0) {
			optimal_partitions.push_back(partition);
		} else {
			// get size of partition one level below, i.e. number of nodes in original graph
			unsigned int original_number_of_nodes =
					optimal_partitions[hierarchy - 1].element_count();
			// create new empty partition of this size
			P partition_original_nodes(original_number_of_nodes);

			// loop over nodes one level below
			int old_comm, new_comm;
			for (unsigned int id = 0; id < original_number_of_nodes; id++) {
				// get the communities for each node one level below
				old_comm = optimal_partitions[hierarchy - 1].find_set(id);
				// use this as node_id in the current partition as old community id
				//is equivalent to new node id and read out new community id
				new_comm = partition.find_set(old_comm);
				// include pair (node, new community) id in the newly created partition
				partition_original_nodes.add_node_to_set(id, new_comm);
			}
			optimal_partitions.push_back(partition_original_nodes);
		}

		// get number of communities
		int num_comm = partition.set_count();

		// Create graph from partition
		T reduced_graph;
		for (int i = 0; i < num_comm; i++) {
			reduced_graph.addNode();
		}

		//Need a map set_id > node_in_reduced_graph_
		W reduced_weight_map(reduced_graph);

		// Find between community total weights by checking
		// Edges within graph
		for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
			int comm_of_node_u = partition.find_set(graph.id(graph.u(edge)));
			int comm_of_node_v = partition.find_set(graph.id(graph.v(edge)));

			float weight = weights[edge];

			Node node_of_comm_u = reduced_graph.nodeFromId(comm_of_node_u);
			Node node_of_comm_v = reduced_graph.nodeFromId(comm_of_node_v);

			Edge edge_in_reduced_graph = lemon::findEdge(reduced_graph,
					node_of_comm_u, node_of_comm_v);

			if (edge_in_reduced_graph == lemon::INVALID) {
				edge_in_reduced_graph = reduced_graph.addEdge(node_of_comm_u,
						node_of_comm_v);
				reduced_weight_map[edge_in_reduced_graph] = weight;
			} else {
				reduced_weight_map[edge_in_reduced_graph] += weight;
			}
		}

		return cliques::find_optimal_partition_louvain_with_gain<P>(
				reduced_graph, reduced_weight_map, compute_quality,
				compute_quality_diff, optimal_partitions);
	}

	return current_quality;
}

}

#endif
