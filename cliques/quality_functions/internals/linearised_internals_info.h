#pragma once

#include <cliques/helpers/graphhelpers.h>

namespace clq {
// Define internal structure to carry statistics for combinatorial stability

struct LinearisedInternalsInfo {

	// typedef for convenience
	typedef std::vector<double> range_map;

	unsigned int num_nodes; // number of nodes in graph
	range_map comm_w_in; // weight inside each community

	std::vector<double> node_weight_to_communities; //mapping: node weight to each community
	std::vector<unsigned int> neighbouring_communities_list; // associated list of neighbouring communities

	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// simple constructor, no partition given; graph should be the original graph in this case
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	template<typename G, typename M>
	LinearisedInternalsInfo(G &graph, M &weights) :
		num_nodes(lemon::countNodes(graph)), comm_w_in(num_nodes, 0),
				node_weight_to_communities(num_nodes, 0),
				neighbouring_communities_list() {
		for (unsigned int i = 0; i < num_nodes; ++i) {
			typename G::Node temp_node = graph.nodeFromId(i);
			comm_w_in[i] = find_weight_selfloops(graph, weights, temp_node);
		}
	}
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// full constructor with reference to partition
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	template<typename G, typename M, typename P>
	LinearisedInternalsInfo(G &graph, M &weights, P &partition) :
		num_nodes(lemon::countNodes(graph)), comm_w_in(num_nodes, 0),
				node_weight_to_communities(lemon::countNodes(graph), 0),
				neighbouring_communities_list() {

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

			// in case the weight stems from within the community add to internal weights
			if (comm_of_node_u == comm_of_node_v) {
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
		LinearisedInternalsInfo &internals, P &partition) {
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
		if (graph.u(e) != graph.v(e)) {
			double edge_weight = weights[e];
			typename G::Node opposite_node = graph.oppositeNode(node, e);
			int comm_node = partition.find_set(graph.id(opposite_node));
			if (internals.node_weight_to_communities[comm_node] == 0) {
				internals.neighbouring_communities_list.push_back(comm_node);
			}
			internals.node_weight_to_communities[comm_node] += edge_weight;
		}
	}

	internals.comm_w_in[comm_id] -= 2
			* internals.node_weight_to_communities[comm_id]
			+ find_weight_selfloops(graph, weights, node);
	//clq::output("in", internals.comm_w_in[comm_id]);

	partition.isolate_node(node_id);
}

/**
 @brief  insert a node into the best set & update internals
 */
template<typename G, typename M, typename P>
void insert_and_update_internals(G &graph, M &weights, typename G::Node node,
		LinearisedInternalsInfo &internals, P &partition, int best_comm) {
	int node_id = graph.id(node);
	// insert node to partition/bookkeeping

	internals.comm_w_in[best_comm] += 2
			* internals.node_weight_to_communities[best_comm]
			+ find_weight_selfloops(graph, weights, node);
	partition.add_node_to_set(node_id, best_comm);
//	    clq::output("in", internals.comm_w_in[best_comm]);
//	    clq::print_partition_line(partition);
}

}
