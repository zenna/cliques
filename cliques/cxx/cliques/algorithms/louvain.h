/* Copyright (c) Z Tavares, M Schaub - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <limits>
#include <vector>
#include <iostream>
#include <ctime>

#include <lemon/concepts/graph.h>

#include <cliques/helpers.h>
#include <cliques/graphhelpers.h>
#include <cliques/algorithms/internals/internals.h>

// TODO: Louvain gets noticably slower on second iteration
// TODO: It Doesn't even when map of different graph is passed, this is wrong
// TODO: On profiling, isolate_and_update_internals and find_selfloops seems to be bottleneck

namespace cliques {

/**
 @brief  Louvain method - greedy algorithm to find community structure of a network.

 The second phase of the algorithm consists in building a new network
 whose nodes are now the communities found during the first phase.
 To do so, the weights of the links between the new nodes are given by
 the sum of the weight of the links between nodes in the corresponding
 two communities.
 Links between nodes of the same community lead to self-loops for this
 community in the new network. Once this second phase is completed,
 it is then possible to reapply the first phase of the algorithm
 to the resulting weighted network and to iterate.

 @param[in]  graph     graph to partition
 @param[in]  weights   EdgeMap of edge weights
 @param[in]  compute_quality     partition quality functor
 @param[out] optimal_partitions     optimal partitions found, last in vector is overall best
 */
template<typename P, typename T, typename W, typename QF, typename QFDIFF, typename Logger>
double find_optimal_partition_louvain_with_gain(T &graph, W &weights,
		QF compute_quality, QFDIFF compute_quality_diff,
		std::vector<P> &optimal_partitions, Logger &log) {

	typedef typename T::Node Node;
	typedef typename T::Edge Edge;
	typedef typename T::NodeIt NodeIt;
	typedef typename T::EdgeIt EdgeIt;
	typedef typename T::IncEdgeIt IncEdgeIt;

	auto internals = cliques::gen(compute_quality, graph, weights);
	P partition(lemon::countNodes(graph));
	partition.initialise_as_singletons();
	double minimum_improve = 0.000000001;
	double current_quality = compute_quality(internals);
	bool one_level_end = false;
	double old_quality = current_quality;
	bool did_nodes_move = false;

	// Randomise the looping over nodes
	// initialise random number generator
	srand(std::time(0));
	// create vector to shuffle
	std::vector<Node> nodes_ordered_randomly;
	for (NodeIt temp_node(graph); temp_node != lemon::INVALID; ++temp_node) {
		nodes_ordered_randomly.push_back(temp_node);
	}
	// Randomly shuffle nodes
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

			// HIJACK ---------------------------------------
			int hierarchy = optimal_partitions.size();
			if (hierarchy == 0) {
				log.log(partition);
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
				log.log(partition_original_nodes);
			}
			// ----------------------------------------------

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

	int hierarchy = optimal_partitions.size();
	if (one_level_end == true) {
		// Compile P into original partition size
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
		int i = 0;
		for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
			int comm_of_node_u = partition.find_set(graph.id(graph.u(edge)));
			int comm_of_node_v = partition.find_set(graph.id(graph.v(edge)));
			i++;

			float weight = weights[edge];

			Node node_of_comm_u = reduced_graph.nodeFromId(comm_of_node_u);
			Node node_of_comm_v = reduced_graph.nodeFromId(comm_of_node_v);

			// TODO: findEdge is slow!
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
				compute_quality_diff, optimal_partitions, log);
	}
	if (hierarchy == 0) {
		optimal_partitions.push_back(partition);
	}
	return current_quality;
}

}
