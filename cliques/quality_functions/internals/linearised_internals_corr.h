#pragma once

#include <cliques/helpers/graphhelpers.h>

namespace clq {
// Define internal structure to carry statistics for correlation based normalised stability

struct LinearisedInternalsCorr {

	// typedef for convenience
	typedef std::vector<double> range_map;

	unsigned int num_nodes; // number of nodes in graph
	unsigned int num_nodes_init; // number of nodes in original graph
	double two_m;
	range_map null_model;
	range_map node_to_w; // weighted loss of each node
	range_map comm_loss; // loss for each community (second matrix / null model term)
	range_map comm_w_in; // gain for each community (first matrix / gain term)

	std::vector<double> node_weight_to_communities; //mapping: node weight to each community
	std::vector<unsigned int> neighbouring_communities_list; // associated list of neighbouring communities

	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// simple constructor, no partition given; graph should be the original graph in this case
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	template<typename G, typename M>
	LinearisedInternalsCorr(G &graph, M &weights,
			std::vector<double> null_model_vec) :
		num_nodes(lemon::countNodes(graph)), null_model(num_nodes, 0),
				node_to_w(num_nodes, 0), comm_loss(num_nodes, 0), comm_w_in(
						num_nodes, 0),
				node_weight_to_communities(num_nodes, 0),
				neighbouring_communities_list() {

		num_nodes_init = num_nodes; // get original number of nodes

		for (unsigned int i = 0; i < num_nodes_init; ++i) {
			typename G::Node temp_node = graph.nodeFromId(i);
			node_to_w[i] = null_model_vec[i] / sqrt(null_model_vec[i] * (1
					- null_model_vec[i]));
			comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node);
			comm_loss[i] = null_model_vec[i] / sqrt(null_model_vec[i] * (1
					- null_model_vec[i]));
		}
	}
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// full constructor with reference to partition
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	template<typename G, typename M, typename P>
	LinearisedInternalsCorr(G &graph, M &weights, P &partition,
			P &partition_init, std::vector<double> null_model_vec) :
		num_nodes(lemon::countNodes(graph)), null_model(null_model_vec),
				node_to_w(num_nodes, 0), comm_loss(num_nodes, 0), comm_w_in(
						num_nodes, 0),
				node_weight_to_communities(num_nodes, 0),
				neighbouring_communities_list() {

		typedef typename G::EdgeIt EdgeIt;
		num_nodes_init = partition_init.element_count(); // get original number of nodes

		for (unsigned int i = 0; i < num_nodes_init; ++i) {
			int old_comm_id = partition_init.find_set(i); // get community of original node
			// old_comm_id is equal to node id in new graph
			node_to_w[old_comm_id] += null_model[i] / sqrt(
					null_model_vec[i] * (1 - null_model_vec[i]));
			int new_comm_id = partition.find_set(old_comm_id);
			// compute loss based on original graph
			comm_loss[new_comm_id] += null_model[i] / sqrt(
					null_model_vec[i] * (1 - null_model_vec[i]));
		}

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

			// if selfloop, only half of the weight has to be considered
			if (node_u_id == node_v_id) {
				weight = weight / 2;
			}

			if (comm_of_node_u == comm_of_node_v) {
				// in case the weight stems from within the community add to internal weights
				comm_w_in[comm_of_node_u] += 2 * weight;
			}
		}

	}
};

/**
 @brief  isolate a node into its singleton set & update internals
 */
template<typename G, typename M, typename P>
void isolate_and_update_internals(G &graph, M &weights, typename G::Node node,
		LinearisedInternalsCorr &internals, P &partition) {
	int node_id = graph.id(node);
	int comm_id = partition.find_set(node_id);

	// reset weights
	while (!internals.neighbouring_communities_list.empty()) {
		//clq::output("empty");
		unsigned int old_neighbour =
				internals.neighbouring_communities_list.back();
		internals.neighbouring_communities_list.pop_back();
		internals.node_weight_to_communities[old_neighbour] = 0;
	}

	// get weights from node to each community
	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		// check that you do not get a self-loop
		if (graph.u(e) != graph.v(e)) {
			// get the edge weight
			double edge_weight = weights[e];
			// get the other node
			typename G::Node opposite_node = graph.oppositeNode(node, e);
			// get community id of the other node
			int comm_node = partition.find_set(graph.id(opposite_node));
			// check if we have seen this community already
			if (internals.node_weight_to_communities[comm_node] == 0) {
				internals.neighbouring_communities_list.push_back(comm_node);
			}
			// add weights to vector
			internals.node_weight_to_communities[comm_node] += edge_weight;
		}
	}
//	clq::print_collection(internals.node_weight_to_communities);
//	clq::print_partition_line(partition);
//	clq::output("loss", internals.comm_loss[comm_id]);
	internals.comm_loss[comm_id] -= internals.node_to_w[node_id];
//	clq::output("loss", internals.comm_loss[comm_id]);
//	clq::output("in", internals.comm_w_in[comm_id]);
	internals.comm_w_in[comm_id] -= 2
			* internals.node_weight_to_communities[comm_id]
			+ find_weight_selfloops(graph, weights, node);
//	clq::output("in", internals.comm_w_in[comm_id]);

	partition.isolate_node(node_id);
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename G, typename M, typename P>
void insert_and_update_internals(G &graph, M &weights, typename G::Node node,
		LinearisedInternalsCorr &internals, P &partition, int best_comm) {
	// node id and std dev
	int node_id = graph.id(node);

	// update loss
	internals.comm_loss[best_comm] += internals.node_to_w[node_id];
	// update gain
	internals.comm_w_in[best_comm] += 2
			* internals.node_weight_to_communities[best_comm]
			+ find_weight_selfloops(graph, weights, node);

	partition.add_node_to_set(node_id, best_comm);
	//    clq::output("in", internals.comm_w_in[best_comm], "tot",internals.comm_w_tot[best_comm], "nodew", internals.node_to_w[best_comm], "2m", internals.two_m);
	//    clq::print_partition_line(partition);
}

}
