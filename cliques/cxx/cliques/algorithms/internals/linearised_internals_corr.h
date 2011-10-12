#pragma once

#include <cliques/graphhelpers.h>

namespace cliques {
// Define internal structure to carry statistics for correlation based normalised stability

struct LinearisedInternalsCorr {

	// typedef for convenience
	typedef std::vector<double> range_map;

	unsigned int num_nodes; // number of nodes in graph
	double two_m; // 2 times total weight
	range_map node_to_w; // weighted degree of each node
	range_map comm_loss; // loss for each community (second/null model term)
	range_map comm_w_in; // gain for each community (first/gain term)

	std::vector<double> node_weight_to_communities; //mapping: node weight to each community
	std::vector<unsigned int> neighbouring_communities_list; // associated list of neighbouring communities

	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// simple constructor, no partition given; graph should be the original graph in this case
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	template<typename G, typename M>
	LinearisedInternalsCorr(G &graph, M &weights) :
		num_nodes(lemon::countNodes(graph)), comm_loss(num_nodes, 0),
				comm_w_in(num_nodes, 0), node_to_w(num_nodes, 0),
				node_weight_to_communities(num_nodes, 0),
				neighbouring_communities_list() {

		two_m = 2 * find_total_weight(graph, weights); // total graph weight

		for (unsigned int i = 0; i < num_nodes; ++i) {
			typename G::Node temp_node = graph.nodeFromId(i);
			node_to_w[i] = find_weighted_degree(graph, weights, temp_node);
			// correlation => normalised by variance (time 0)
			double var_times_2m = node_to_w[i] * (1 - node_to_w[i] / two_m);
			comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node)
					/ var_times_2m;
			comm_loss[i] = node_to_w[i] / sqrt(node_to_w[i] * (two_m
					- node_to_w[i]));
		}
	}
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// full constructor with reference to partition
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	template<typename G, typename M, typename P>
	LinearisedInternalsCorr(G &graph, M &weights, P &partition) :
		num_nodes(lemon::countNodes(graph)), node_to_w(num_nodes, 0),
				comm_loss(num_nodes, 0), comm_w_in(num_nodes, 0),
				node_weight_to_communities(num_nodes, 0),
				neighbouring_communities_list() {

		two_m = 2 * find_total_weight(graph, weights);

		for (unsigned int i = 0; i < num_nodes; ++i) {
			typename G::Node temp_node = graph.nodeFromId(i);
			node_to_w[i] = find_weighted_degree(graph, weights, temp_node);
		}

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

			// if selfloop, only half of the weight has to be considered
			if (node_u_id == node_v_id) {
				weight = weight / 2;
			}

			double sigma_u_times_2m = sqrt(node_to_w[node_u_id] * (two_m
					- node_to_w[node_u_id]));
			double sigma_v_times_2m = sqrt(node_to_w[node_v_id] * (two_m
					- node_to_w[node_v_id]));

			// add weight to total weight of community
			comm_loss[comm_of_node_u] += weight / sigma_u_times_2m;
			comm_loss[comm_of_node_v] += weight / sigma_v_times_2m;

			if (comm_of_node_u == comm_of_node_v) {
				// in case the weight stems from within the community add to internal weights
				comm_w_in[comm_of_node_u] += 2 * weight * two_m
						/ (sigma_u_times_2m * sigma_v_times_2m);
			}
		}

	}
};

/**
 @brief  isolate a node into its singleton set & update internals
 */
template<typename G, typename M, typename I, typename P>
void isolate_and_update_internals(G &graph, M &weights, typename G::Node node,
		I &internals, P &partition) {
	int node_id = graph.id(node);
	int comm_id = partition.find_set(node_id);

	// reset weights
	while (!internals.neighbouring_communities_list.empty()) {
		//cliques::output("empty");
		unsigned int old_neighbour =
				internals.neighbouring_communities_list.back();
		internals.neighbouring_communities_list.pop_back();
		internals.node_weight_to_communities[old_neighbour] = 0;
	}

	// std deviation of node itself
	double sigma_i_times_2m = sqrt(internals.node_to_w[node_id]
			* (internals.two_m - internals.node_to_w[node_id]));

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
			int other_node_id = graph.id(opposite_node);
			double sigma_j_times_2m = sqrt(internals.node_to_w[other_node_id]
					* (internals.two_m - internals.node_to_w[other_node_id]));
			internals.node_weight_to_communities[comm_node] += edge_weight
					* internals.two_m / (sigma_i_times_2m * sigma_j_times_2m);
		}
	}
	//cliques::print_collection(internals.node_weight_to_communities);
	internals.comm_loss[comm_id] -= internals.node_to_w[node_id]
			/ sigma_i_times_2m;
	//cliques::output("in", internals.comm_w_in[comm_id]);
	internals.comm_w_in[comm_id] -= 2
			* internals.node_weight_to_communities[comm_id]
			+ find_weight_selfloops(graph, weights, node) * internals.two_m
					/ (sigma_i_times_2m * sigma_i_times_2m);
	//cliques::output("in", internals.comm_w_in[comm_id]);

	partition.isolate_node(node_id);
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename G, typename M, typename I, typename P>
void insert_and_update_internals(G &graph, M &weights, typename G::Node node,
		I &internals, P &partition, int best_comm) {
	// node id and std dev
	int node_id = graph.id(node);
	double sigma_i_times_2m = sqrt(internals.node_to_w[node_id]
			* (internals.two_m - internals.node_to_w[node_id]));

	// update loss
	internals.comm_loss[best_comm] += internals.node_to_w[node_id]
			/ sigma_i_times_2m;
	// update gain
	internals.comm_w_in[best_comm] += 2
			* internals.node_weight_to_communities[best_comm]
			+ find_weight_selfloops(graph, weights, node) * internals.two_m
					/ (sigma_i_times_2m * sigma_i_times_2m);

	partition.add_node_to_set(node_id, best_comm);
	//    cliques::output("in", internals.comm_w_in[best_comm], "tot",internals.comm_w_tot[best_comm], "nodew", internals.node_to_w[best_comm], "2m", internals.two_m);
	//    cliques::print_partition_line(partition);
}

}
