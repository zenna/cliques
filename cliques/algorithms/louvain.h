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
#include <cliques/algorithms/internals/generators.h>


// TODO: Separate out louvain into smaller functions
// TODO: update random shuffling to not depend on time random seed

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
template<typename P, typename T, typename W, typename QF, typename QFDIFF,
		typename Logger>
double find_optimal_partition_louvain(T &graph, W &weights, std::vector<
		double> null_model_vec, QF compute_quality,
		QFDIFF compute_quality_diff, P initial_partition,
		std::vector<P> &optimal_partitions, double minimum_improve, Logger log,
		bool log_switch = false) {

	typedef typename T::Node Node;
	typedef typename T::Edge Edge;
	typedef typename T::NodeIt NodeIt;
	typedef typename T::EdgeIt EdgeIt;
	typedef typename T::IncEdgeIt IncEdgeIt;

	P partition(lemon::countNodes(graph));
	partition = initial_partition;
	P partition_init = initial_partition;
	if (!optimal_partitions.empty()) {
		partition_init = optimal_partitions.back();
	} auto internals = cliques::gen(compute_quality, graph, weights, partition, partition_init, null_model_vec);

	//double minimum_improve = 0.000000001; //1e-9
	double current_quality = compute_quality(internals);
	//cliques::output("current_quality", current_quality);
	bool one_level_end = false;
	double old_quality = current_quality;
	bool did_nodes_move = false;

	// Randomise the looping over nodes
	// initialise random number generator
	//srand(std::time(0));
	//TODO find neat way of doing this once only...
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

			unsigned int num_neighbour_communities =
					internals.neighbouring_communities_list.size();

			// loop over neighbouring communities
			for (unsigned int k = 0; k < num_neighbour_communities; ++k) {

				unsigned int comm_id_neighbour =
						internals.neighbouring_communities_list[k];
				double gain = compute_quality_diff(internals,
						comm_id_neighbour, node_id);

				if (gain > best_gain) {
					best_comm = comm_id_neighbour;
					best_gain = gain;
					// avoid not necessary movements, place node in old community if possible
				} else if (gain == best_gain && comm_id == comm_id_neighbour) {
					best_comm = comm_id;
				}

			}

			insert_and_update_internals(graph, weights, n1, internals,
					partition, best_comm);

			// Logging ---------------------------------------
			if (log_switch) {
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
						old_comm = optimal_partitions[hierarchy - 1].find_set(
								id);
						// use this as node_id in the current partition as old community id
						//is equivalent to new node id and read out new community id
						new_comm = partition.find_set(old_comm);
						// include pair (node, new community) id in the newly created partition
						partition_original_nodes.add_node_to_set(id, new_comm);
					}
					log.log(partition_original_nodes);
				}
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
	std::map<int, int> new_comm_id_to_old_comm_id = partition.normalise_ids();

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

		// Create graph from partition
		T reduced_graph;
		W reduced_weight_map(reduced_graph);
		create_reduced_graph_from_partition(reduced_graph, reduced_weight_map,
				graph, weights, partition, new_comm_id_to_old_comm_id,
				internals);
		if(log_switch){
			std::string filename("intermediate_graphs_");
			std::stringstream hierarchy_level_sstream;
			hierarchy_level_sstream << hierarchy + 1;
			std::string hierarchy_level(hierarchy_level_sstream.str());
			filename += hierarchy_level_sstream.str();
			write_edgelist_weighted(filename, reduced_graph, reduced_weight_map);
		}
		P singleton_partition(lemon::countNodes(reduced_graph));
		singleton_partition.initialise_as_singletons();

		return cliques::find_optimal_partition_louvain<P>(reduced_graph,
				reduced_weight_map, null_model_vec, compute_quality,
				compute_quality_diff, singleton_partition, optimal_partitions,
				minimum_improve, log, log_switch);
	}
	if (hierarchy == 0) {
		optimal_partitions.push_back(partition);
	}
	return current_quality;
}

}
