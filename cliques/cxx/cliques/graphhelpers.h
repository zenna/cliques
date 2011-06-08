/* Copyright (c) Zenna Tavares & M. Schaub - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <lemon/maps.h>
#include <lemon/core.h>

namespace cliques {

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

template <typename G>
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
		// TODO make aware that this is counting self loops twice..
		degree += weights[e];
	}
	return degree;
}

template<typename G, typename M, typename NO>
double find_weight_selfloops(G &graph, M &weights, NO node) {
	typename G::Edge edge = lemon::findEdge(graph, node, node);
	if (edge == lemon::INVALID) {
		return 0.0;
	}
	else {
		return 2*double(weights[edge]);
	}
}

template<typename G, typename M>
float find_total_weight(G &graph, M &weights) {
	double total_weight = 0.0;
	for (typename G::EdgeIt e(graph); e != lemon::INVALID; ++e) {
		total_weight = total_weight + weights[e];
	}
	return total_weight;
}

template<typename G, typename P, typename W, typename NO>
std::map<int, double> find_weight_node_to_communities(G &graph, P &partition,
		W &weights, NO node) {
	std::map<int, double> community_to_weight;

	for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
		double edge_weight = weights[e];
		NO opposite_node = graph.oppositeNode(node, e);
		int comm_node = partition.find_set(graph.id(opposite_node));
		community_to_weight[comm_node] += edge_weight;
	}

	return community_to_weight;
}

template<typename G>
void read_edgelist(G &graph, std::string filename) {
	std::ifstream maxima_file(filename.c_str());
	std::string line;
	std::string maxima;

	if (!maxima_file.is_open()) {
		std::cout << "couldn't open file" << std::endl;
		exit(1);
	}

	typedef typename G::Node Node;
	std::map<int, Node> id_to_node;
	while (std::getline(maxima_file, line)) {
		std::stringstream lineStream(line);
		std::set<int> current_maxima;

		std::getline(lineStream, maxima, ' ');
		int node1_id = atoi(maxima.c_str());
		std::getline(lineStream, maxima, ' ');
		int node2_id = atoi(maxima.c_str());

		typename std::map<int, Node>::iterator itr = id_to_node.find(node1_id);
		Node node1, node2;

		if (itr == id_to_node.end()) {
			node1 = graph.addNode();
			id_to_node[node1_id] = node1;
		} else {
			node1 = itr->second;
		}
		itr = id_to_node.find(node2_id);
		if (itr == id_to_node.end()) {
			node2 = graph.addNode();
			id_to_node[node2_id] = node2;
		} else {
			node2 = itr->second;
		}
		graph.addEdge(node1, node2);
	}
	maxima_file.close();
}

template<typename G, typename E>
void read_edgelist_weighted(std::string filename, G &graph, E &weights) {
	// initialise input stream and strings for readout
	std::ifstream maxima_file(filename.c_str());
	std::string line;
	std::string maxima;

	// check if file is open
	if (!maxima_file.is_open()) {
		std::cout << "couldn't open file" << std::endl;
		exit(1);
	}

	// define Node class for convenience
	typedef typename G::Node Node;
	// mapping from id to node
	std::map<int, Node> id_to_node;

	//readout contents from maxima_file into string, line by line
	while (std::getline(maxima_file, line)) {

		std::stringstream lineStream(line);
		//readout node id and weights
		std::getline(lineStream, maxima, ' ');
		int node1_id = atoi(maxima.c_str());
		std::getline(lineStream, maxima, ' ');
		int node2_id = atoi(maxima.c_str());
		std::getline(lineStream, maxima, ' ');
		float weight = atof(maxima.c_str());

		typename std::map<int, Node>::iterator itr = id_to_node.find(node1_id);
		Node node1, node2;

		// If the node is not in the map
		// then create node and add to map
		if (itr == id_to_node.end()) {
			node1 = graph.addNode();
			id_to_node[node1_id] = node1;
		} else {
			node1 = itr->second;
		}

		// same for node 2
		itr = id_to_node.find(node2_id);
		if (itr == id_to_node.end()) {
			node2 = graph.addNode();
			id_to_node[node2_id] = node2;
		} else {
			node2 = itr->second;
		}

		//std::cout << "adding " << graph.id(node1) << " - " << graph.id(node2) << std::endl;
		typename G::Edge edge = graph.addEdge(node1, node2);
		weights.set(edge, weight);
	}

	maxima_file.close();
}

/**
 @brief  Checks whether a partition is still connected

 A connected partition is one whereby all nodes within any one group
 must be mutually accessible to all others within the group without
 having to traverse nodes outside that group.

 @param[in]  partition       partition to evaluate
 */
template <typename P>
bool is_partition_connected(P &partition) {
	return true;
}

}
