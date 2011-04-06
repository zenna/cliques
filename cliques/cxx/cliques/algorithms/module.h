/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_MODULE_H
#define CLIQUES_MODULE_H

#include <limits>
#include <vector>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/structures/partition.h>
#include <cliques/structures/common.h>

namespace cliques {

/**
 @brief  Louvain method - greedy algorithm to find community structure of a network.

 This is a fast algorithm to find communities within a network
 ref: Fast Unfolding Of Communities in large networks, Blondel et al. 2008

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

/* Implementation for Louvain algorithm for "modularity-like" quality functions, where the change
 * in quality can be evaluated efficiently.
 */

template<typename P, typename T, typename W, typename QF, typename QFDIFF>
void find_optimal_partition_louvain_with_gain(T &graph, W &weights,
		QF quality_function, QFDIFF quality_function_diff,
		std::vector<P> &optimal_partitions) {

	// Start: Create singleton partition from graph
	P partition;
	for (typename T::NodeIt n(graph); n != lemon::INVALID; ++n) {
		partition.add_element(graph.id(n));
	}

	// some bookkeeping	to keep matters simple
	// minimum improvement
	double minimum_improve = 0.000001; //0;
	// total weight of the links
	double two_m = find_total_weight(graph, weights);
	// number_of_nodes
	unsigned int number_of_nodes = lemon::countNodes(graph);
	// mapping form node to node weight
	lemon::RangeMap<double> node_to_w =
			lemon::RangeMap<double>(number_of_nodes);
	// mapping from node number to community number
	lemon::RangeMap<unsigned int> node_to_comm = lemon::RangeMap<unsigned int>(
			number_of_nodes, 0);
	// mapping #community to total weight of community
	lemon::RangeMap<double> comm_w_tot = lemon::RangeMap<double>(
			number_of_nodes, 0);
	// mapping #community to weight within community
	lemon::RangeMap<double> comm_w_in = lemon::RangeMap<double>(
			number_of_nodes, 0);

	// initialise bookkeeping variables..
	for (int i = 0; i < number_of_nodes; i++) {
		typename T::Node temp_node = graph.nodeFromId(i);
		comm_w_tot[i] = node_to_w[i] = find_weighted_degree(graph, weights,
				temp_node);
		node_to_comm[i] = i;
		comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node);
	}

	double current_quality = quality_function(comm_w_tot, comm_w_in, two_m);

	//TODO: Randomise the looping over nodes

	bool one_level_end = false;
	unsigned int number_of_moves = 0;
	double old_quality = current_quality;
	//one level louvain
	do {
		// initialise number of moves, get current quality
		number_of_moves = 0;
		old_quality = current_quality;

		// Loop over all nodes
		for (typename T::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {

			// get node id and comm id
			unsigned int node_id = graph.id(n1);
			// TODO this has to be done with real partition
			unsigned int comm_id = node_to_comm[node_id];

			// get weights from node to each community
			lemon::RangeMap<double> node_weight_to_communities =
					cliques::find_weight_node_to_communities(graph, partition,
							weights);

			//TODO: make inline function for this
			// -----------------------------------------
			// remove node from partition/bookkeeping
			comm_w_tot[comm_id] -= node_to_w[node_id];
			comm_w_in[comm_id] -= 2 * node_weight_to_communities[comm_id]
					+ find_weight_selfloops(graph, weights, n1);
			// TODO this has to be done with real partition
			node_to_comm[node_id] = -1;
			// ------------------------------------------


			//default option for re-inclusion of node
			unsigned int best_comm = comm_id;
			double best_gain = 0;

			// TODO check if better to loop over all neighbouring  communities,
			// find all neighbouring communities

			// loop over all neighbouring nodes
			for (typename T::IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
				typename T::Node n2 = graph.oppositeNode(n1, e);
				if (n1 != n2) {
					// get neighbour node id and neighbour community id
					unsigned int node_id_neighbour = graph.id(n2);
					unsigned int comm_id_neighbour =
							node_to_comm[node_id_neighbour];
					double gain = quality_function_diff(
							comm_w_tot[comm_id_neighbour],
							node_weight_to_communities[comm_id_neighbour],
							two_m, node_to_w[node_id]);
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

			//TODO: make inline function for this
			// -----------------------------------------
			// insert node to partition/bookkeeping
			comm_w_tot[best_comm] += node_to_w[node_id];
			comm_w_in[best_comm] += 2 * node_weight_to_communities[best_comm]
					+ find_weight_selfloops(graph, weights, n1);
			// TODO this has to be done with real partition
			node_to_comm[node_id] = best_comm;
			// ------------------------------------------
			if (best_comm != comm_id) {
				number_of_moves++;
			}

		}

		if (number_of_moves > 0) {
			// compute new quality
			current_quality = quality_function(comm_w_tot, comm_w_in, two_m);
			one_level_end = true;
		}
		// If there was a movement of the nodes AND quality increases make another run
	} while (number_of_moves > 0 && (current_quality - old_quality)
			> minimum_improve);

	// temp. TODO change this
	// create partition from node_to_comm_vector

	// The second phase of the algorithm consists in building a new network
	// whose nodes are now the communities found during the first phase.
	// To do so, the weights of the links between the new nodes are given by
	// the sum of the weight of the links between nodes in the corresponding
	// two communities.
	// Links between nodes of the same community lead to self-loops for this
	// community in the new network. Once this second phase is completed,
	// it is then possible to reapply the first phase of the algorithm
	// to the resulting weighted network and to iterate.
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
			P partition_original_nodes;// = P(original_number_of_nodes);

			// loop over nodes one level below
			int old_comm, new_comm;
			for (int id = 0; id < original_number_of_nodes; id++) {
				// get the communities for each node one level below
				old_comm = optimal_partitions[hierarchy - 1].find_set(id);
				// use this as node_id in the current partition as old community id
				//is equivalent to new node id and read out new community id
				new_comm = partition.find_set(old_comm);
				// include pair (node, new community) id in the newly created partition
				//partition_original_nodes
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
		W reduced_weight_map(reduced_graph);

		//Need a map set_id > node_in_reduced_graph_

		// Find between community total weights by checking
		// Edges within graph
		for (typename T::EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
			int comm_of_node_u = partition.find_set(graph.id(graph.u(edge)));
			int comm_of_node_v = partition.find_set(graph.id(graph.v(edge)));

			float weight = weights[edge];

			typename T::Node node_of_comm_u = reduced_graph.nodeFromId(comm_of_node_u);
			typename T::Node node_of_comm_v = reduced_graph.nodeFromId(comm_of_node_v);

			typename T::Edge edge_in_reduced_graph = lemon::findEdge(
					reduced_graph, comm_of_node_u, comm_of_node_v);

			if (edge_in_reduced_graph == lemon::INVALID) {
				reduced_graph.addEdge(node_of_comm_u);
				reduced_graph.addEdge(node_of_comm_v);
			}
			reduced_weight_map[edge_in_reduced_graph] += weight;
		}

		cliques::find_optimal_partition_louvain<P>(graph, reduced_weight_map,
				quality_function, quality_function_diff, optimal_partitions);

	}

	// Compile P into original partition size
	short hierarchy = optimal_partitions.size();
	if (hierarchy == 0) {
		optimal_partitions.push_back(partition);
	} else {
		// get size of partition one level below, i.e. number of nodes in original graph
		unsigned int original_number_of_nodes = optimal_partitions[hierarchy
				- 1].element_count();
		// create new empty partition of this size
		P partition_original_nodes;// = P(original_number_of_nodes);

		// loop over nodes one level below
		int old_comm, new_comm;
		for (int id = 0; id < original_number_of_nodes; id++) {
			// get the communities for each node one level below
			old_comm = optimal_partitions[hierarchy - 1].find_set(id);
			// use this as node_id in the current partition as old community id
			//is equivalent to new node id and read out new community id
			new_comm = partition.find_set(old_comm);
			// include pair (node, new community) id in the newly created partition
			//partition_original_nodes
		}
		optimal_partitions.push_back(partition_original_nodes);
	}

}

}// end namespace cliques

#endif
