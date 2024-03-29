/* Copyright (c) Zenna Tavares & M. Schaub - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <lemon/maps.h>
#include <lemon/core.h>
#include <math.h>

namespace clq {

template<typename G>
float A(G &graph, int node1_id, int node2_id) {
	typename G::Node n1 = graph.nodeFromId(node1_id);
	typename G::Node n2 = graph.nodeFromId(node2_id);

	for (typename G::IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
		if (graph.runningNode(e) == n2) {
			return 1.0;
		}
	}
	return 0.0;
}

template<typename G>
int find_unweighted_degree(G &graph, int node_id) {
	int count = 0;
	typename G::Node n = graph.nodeFromId(node_id);

	for (typename G::IncEdgeIt e(graph, n); e != lemon::INVALID; ++e) {
		++count;
	}
	return count;
}

template<typename G, typename M, typename NO>
float find_weighted_degree(G &graph, M &weights, NO node) {
	double degree = 0.0;
	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			degree += weights[e];
		} else {
			degree += weights[e] / 2; // self-loop is iterated twice with IncEdgeIt
		}
	}
	return degree;
}

template<typename G, typename M, typename NO>
double find_weight_selfloops(G &graph, M &weights, NO node) {
	typename G::Edge edge = lemon::findEdge(graph, node, node);
	if (edge == lemon::INVALID) {
		return 0.0;
	} else {
		return double(weights[edge]);
	}
}

template<typename G, typename M>
float find_total_weight(G &graph, M &weights) {
	double total_weight = 0.0;
	for (typename G::EdgeIt e(graph); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			total_weight = total_weight + weights[e];
		} else {
			total_weight = total_weight + weights[e] / 2;
		}
	}
	return total_weight;
}

template<typename G, typename W>
std::vector<double> create_correlation_graph_from_graph(G &graph, W &weights) {

	// some typedefs
	typedef typename G::EdgeIt EdgeIt;
	typedef typename G::Edge Edge;
	typedef typename G::NodeIt NodeIt;
	typedef typename G::Node Node;

	// get number of nodes and preallocate null_model_vec
	unsigned int num_nodes = lemon::countNodes(graph);
	std::vector<double> null_model_vec(num_nodes, 0);
	double two_m = 0;

	// get weights of each node (temp stored in null_model_vec) and total graph weight
	for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
		int node_u_id = graph.id(graph.u(edge));
		int node_v_id = graph.id(graph.v(edge));

		// weight of edge
		double weight = weights[edge];

		// if selfloop, only half of the weight has to be considered
		if (node_u_id == node_v_id) {
			weight = weight / 2;
		}

		// add weight to node weight and sum up total weight
		null_model_vec[node_u_id] += weight;
		null_model_vec[node_v_id] += weight;
		two_m += 2 * weight;

	}

	// normalise link weights by std deviations
	for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
		int node_u_id = graph.id(graph.u(edge));
		int node_v_id = graph.id(graph.v(edge));

		// weight of edge
		double weight = weights[edge];
		// std deviations: remember that null_model_vec stores weighted degree temporarily
		double sigma_u_times_sqrt_two_m = sqrt(null_model_vec[node_u_id] * (1
				- null_model_vec[node_u_id] / two_m));
		double sigma_v_times_sqrt_two_m = sqrt(null_model_vec[node_v_id] * (1
				- null_model_vec[node_v_id] / two_m));
		// renormalize edge weight
		weights[edge] = weight / (sigma_u_times_sqrt_two_m
				* sigma_v_times_sqrt_two_m);

	}

	// compute null model vector
	for (unsigned int node_id = 0; node_id < num_nodes; ++node_id) {
		null_model_vec[node_id] = null_model_vec[node_id] / two_m;
	}
	//clq::print_collection(null_model_vec);
	return null_model_vec;
}

template<typename G, typename P, typename W, typename I>
void create_reduced_graph_from_partition(G &reduced_graph,
		W &reduced_weight_map, G &graph, W &weights, P &partition, std::map<
				int, int> new_comm_id_to_old_comm_id, I &internals) {

	typedef typename G::EdgeIt EdgeIt;
	typedef typename G::Edge Edge;
	typedef typename G::NodeIt NodeIt;
	typedef typename G::Node Node;

	// get number of communities
	int num_comm = partition.set_count();

	reduced_graph.reserveNode(num_comm);

	// add self-loops in new_graph
	for (int i = 0; i < num_comm; i++) {
		Node comm_node = reduced_graph.addNode();
		Edge e = reduced_graph.addEdge(comm_node, comm_node);
		int old_comm_id = new_comm_id_to_old_comm_id[i];
		reduced_weight_map[e] = internals.comm_w_in[old_comm_id];
	}

	// Find between community weights by checking
	// Edges within graph
	for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
		int comm_of_node_u = partition.find_set(graph.id(graph.u(edge)));
		int comm_of_node_v = partition.find_set(graph.id(graph.v(edge)));

		// internal weights already accounted for
		if (comm_of_node_u == comm_of_node_v) {
			continue;
		}

		double weight = weights[edge];

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
}

// find weight to communities excluding self-loops (to own community)
template<typename G, typename P, typename W, typename NO>
std::map<int, double> find_weight_node_to_communities(G &graph, P &partition,
		W &weights, NO node) {
	std::map<int, double> community_to_weight;

	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			double edge_weight = weights[e];
			NO opposite_node = graph.oppositeNode(node, e);
			int comm_node = partition.find_set(graph.id(opposite_node));
			community_to_weight[comm_node] += edge_weight;
		}
	}
	//clq::print_map(community_to_weight);
	return community_to_weight;
}

// find weight to community excluding self-loops (to own community)
template<typename G, typename P, typename W, typename NO>
double find_weight_node_to_community(G &graph, P &partition, W &weights,
		NO node, int comm_id) {

	double summed_weight = 0.0;
	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		if (graph.u(e) != graph.v(e)) {
			double edge_weight = weights[e];
			NO opposite_node = graph.oppositeNode(node, e);
			int comm_node = partition.find_set(graph.id(opposite_node));
			if (comm_node == comm_id) {
				summed_weight += edge_weight;
			}
		}
	}

	return summed_weight;
}

}
