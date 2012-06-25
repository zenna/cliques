/* Copyright (c) Modified by Zenna Tavares-zennatavares@gmail.com, 2010-2011 */
#pragma once
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <lemon/concepts/graph.h>

namespace clq {

template<typename G, typename M>
void make_fish_graph(G &graph, M &weights, double epsilon, bool many_links) {
	int num_nodes = 8;
	typedef typename G::Node Node;
	typedef typename G::Edge Edge;
	for (int i = 0; i < num_nodes; ++i) {
		graph.addNode();
	}

	// fish head
	for (int i = 0; i < 4; ++i) {
		for (int j = i + 1; j < 4; ++j) {
			Node n1 = graph.nodeFromId(i);
			Node n2 = graph.nodeFromId(i);
			graph.addEdge(n1, n2);
		}
	}

	// tail fin
	graph.addEdge(graph.nodeFromId(4), graph.nodeFromId(5));
	graph.addEdge(graph.nodeFromId(4), graph.nodeFromId(6));
	graph.addEdge(graph.nodeFromId(4), graph.nodeFromId(7));
	graph.addEdge(graph.nodeFromId(5), graph.nodeFromId(6));
	graph.addEdge(graph.nodeFromId(6), graph.nodeFromId(7));

	// now make edge weights
	make_weights_from_edges(graph, weights);

	if (many_links == false) {
		Edge e = graph.addEdge(graph.nodeFromId(3), graph.nodeFromId(4));
		weights[e] = epsilon;
	} else {
		int Ne = 6;
		Edge e1 = graph.addEdge(graph.nodeFromId(3), graph.nodeFromId(4));
		weights[e1] = epsilon / Ne;
		Edge e2 = graph.addEdge(graph.nodeFromId(0), graph.nodeFromId(6));
		weights[e2] = epsilon / Ne;
		Edge e3 = graph.addEdge(graph.nodeFromId(1), graph.nodeFromId(5));
		weights[e3] = epsilon / Ne;
		Edge e4 = graph.addEdge(graph.nodeFromId(2), graph.nodeFromId(7));
		weights[e4] = epsilon / Ne;
		Edge e5 = graph.addEdge(graph.nodeFromId(2), graph.nodeFromId(6));
		weights[e5] = epsilon / Ne;
		Edge e6 = graph.addEdge(graph.nodeFromId(3), graph.nodeFromId(5));
		weights[e6] = epsilon / Ne;
	}

}

template<typename G>
void make_complete_graph(G &graph, int num_nodes) {
	typedef typename G::NodeIt NodeIt;

	for (int i = 0; i < num_nodes; ++i) {
		graph.addNode();
	}

	for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		for (NodeIt n2 = n1; ++n2 != lemon::INVALID;) {
			graph.addEdge(n1, n2);
		}
	}
}

template<typename G, typename M>
void make_complete_graph(G &graph, int num_nodes, M &weights) {
	make_complete_graph(graph, num_nodes);
	make_weights_from_edges(graph, weights);
}

template<typename G>
void make_path_graph(G &graph, int num_nodes) {
	typedef typename G::NodeIt NodeIt;

	for (int i = 0; i < num_nodes; ++i) {
		graph.addNode();
	}
	for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		NodeIt n2 = n1;
		++n2;
		if (n2 == lemon::INVALID) {
			break;
		}
		graph.addEdge(n1, n2);
	}
}

template<typename G, typename M>
void make_path_graph(G &graph, int num_nodes, M &weights) {
	make_path_graph(graph, num_nodes);
	make_weights_from_edges(graph, weights);
}

template<typename G>
void make_ring_graph(G &graph, int num_nodes) {
	typedef typename G::NodeIt NodeIt;

	for (int i = 0; i < num_nodes; ++i) {
		graph.addNode();
	}
	for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
		NodeIt n2 = n1;
		++n2;
		if (n2 == lemon::INVALID) {
			graph.addEdge(n1, NodeIt(graph));
			break;
		}
		graph.addEdge(n1, n2);
	}
}

template<typename G, typename M>
void make_ring_graph(G &graph, int num_nodes, M &weights) {
	make_ring_graph(graph, num_nodes);
	make_weights_from_edges(graph, weights);
}

template<typename G>
void make_hierarchical_benchmark_graph(G &graph, double balance,
		double hierarchy, int num_nodes) {
	for (int i = 0; i < num_nodes; ++i) {
		graph.addNode(i);
	}
}

template<typename G, typename M>
void make_weights_from_edges(G &graph, M &weights) {
	for (typename G::EdgeIt e(graph); e != lemon::INVALID; ++e) {
		weights[e] = 1.0;
	}
}

}