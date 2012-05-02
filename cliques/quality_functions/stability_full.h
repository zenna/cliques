/* Copyright (c) Zenna Tavares, Michael Schaub - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <vector>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/helpers/helpers.h>
#include <cliques/helpers/math.h>

#include <cliques/quality_functions/internals/linearised_internals.h>
#include <cliques/helpers/graphhelpers.h>

namespace cliques {

/**
 @brief  Functor for finding full normalised stability of weighted graph
 */
struct find_full_normalised_stability {

	find_linearised_normalised_stability lin_norm_stability;
	lemon::SmartGraph exp_graph;
	lemon::SmartGraph::EdgeMap<double> exp_graph_weights;
	std::vector<double> minus_t_D_inv_L;
	std::vector<double> node_weighted_degree;
	double m_time;
	double threshold;

	template<typename G, typename M>
	find_full_normalised_stability(G &graph, M &weights, double thres) :
		lin_norm_stability(1), exp_graph(), exp_graph_weights(exp_graph) {

		// get threshold and set time
		threshold = thres;
		m_time =-1;

		typedef typename G::Node Node;
		typedef typename G::NodeIt NodeIt;
		typedef typename G::EdgeIt EdgeIt;
		/////////////////
		// read out graph into matrix
		//////////////////

		// number of nodes N
		int N = lemon::countNodes(graph);
		// -D^-1*L*t == t(B-I)
		minus_t_D_inv_L = std::vector<double>(N * N, 0);
		node_weighted_degree = std::vector<double>(N, 0);

		// get weighted degree of nodes
		for (int i = 0; i < N; ++i) {
			typename G::Node temp_node = graph.nodeFromId(i);
			node_weighted_degree[i] = find_weighted_degree(graph, weights,
					temp_node);
		}

		//initialise matrix and set diagonals to minus identity
		for (int i = 0; i < N * N; ++i) {
			if (i % (N + 1) == 0) {
				minus_t_D_inv_L[i] = -1;
			} else {
				minus_t_D_inv_L[i] = 0;
			}

		}

		// fill in the rest
		for (EdgeIt e(graph); e != lemon::INVALID; ++e) {
			Node u = graph.u(e);
			Node v = graph.v(e);
			int node_id_u = graph.id(u);
			int node_id_v = graph.id(v);
			double weight_uv = weights[e];
			// B_uv
			minus_t_D_inv_L[node_id_v + N * node_id_u] += weight_uv
					/ node_weighted_degree[node_id_v];
			// B_vu
			minus_t_D_inv_L[node_id_u + N * node_id_v] += weight_uv
					/ node_weighted_degree[node_id_u];
		}

		// create graph structure
		//reserve memory space for number of nodes
		exp_graph.reserveNode(N);
		exp_graph.reserveEdge(N + (N * (N - 1)) / 2);

		// add nodes
		for (int i = 0; i < N; ++i) {
			exp_graph.addNode();
		}

	}

	template<typename P>
	double operator ()(P &partition, double markov_time) {

		// create new weight map out of matrix find_stability<lemon::SmartGraph>(partiton,time)
		if (markov_time != m_time) {

			// call expokit
			int N = node_weighted_degree.size();
			std::vector<double> exp_graph_vec = cliques::exp(minus_t_D_inv_L,
					markov_time, N);

			//		cliques::print_collection(exp_graph_vec);


			// reassign markov time
			m_time = markov_time;
			exp_graph.clear();

			// Create new graph
			exp_graph.reserveNode(N);
			exp_graph.reserveEdge(N + (N * (N - 1)) / 2);

			// add nodes
			for (int i = 0; i < N; ++i) {
				exp_graph.addNode();
			}

			// new edges
//			exp_graph_weights = lemon::SmartGraph::EdgeMap<double>(exp_graph);

			for (int i = 0; i < N; ++i) {
				for (int j = i; j < N; ++j) {
					double weight = exp_graph_vec[N * i + j]
							* node_weighted_degree[j];
					if (weight > threshold) {
						lemon::SmartGraph::Edge edge = exp_graph.addEdge(
								exp_graph.nodeFromId(i),
								exp_graph.nodeFromId(j));
						exp_graph_weights.set(edge, weight);
					}
				}
			}
		}

		//		cliques::output("nodes", lemon::countNodes(exp_graph), "edges",lemon::countEdges(exp_graph), "complete", (N * (N-1))/2, N*N, (N*N)/2);
		cliques::LinearisedInternals internals(exp_graph, exp_graph_weights,
				partition);
		return lin_norm_stability(internals);
	}

	template<typename P>
	double operator ()(P &partition, int comm_id, double markov_time) {

		if (markov_time != m_time) {
			// call expokit
			int N = node_weighted_degree.size();

			//        cliques::output("minus");
			//        cliques::print_collection(minus_t_D_inv_L, N);
			std::vector<double> exp_graph_vec = cliques::exp(minus_t_D_inv_L,
					markov_time, N);
			//        cliques::output("exp");
			//        cliques::print_collection(exp_graph_vec, N);

			// create new weight map out of matrix find_stability<lemon::SmartGraph>(partiton,time)

			// reassign markov time
			m_time = markov_time;
			exp_graph.clear();

			// Create new graph
			exp_graph.reserveNode(N);
			exp_graph.reserveEdge(N + (N * (N - 1)) / 2);

			// add nodes
			for (int i = 0; i < N; ++i) {
				exp_graph.addNode();
			}

			// new edges
			for (int i = 0; i < N; ++i) {
				for (int j = i; j < N; ++j) {
					double weight = exp_graph_vec[N * i + j]
							* node_weighted_degree[j];
					if (weight > threshold) {
						lemon::SmartGraph::Edge edge = exp_graph.addEdge(
								exp_graph.nodeFromId(i),
								exp_graph.nodeFromId(j));
						exp_graph_weights.set(edge, weight);
					}
				}
			}
		}

		cliques::LinearisedInternals internals(exp_graph, exp_graph_weights,
				partition);
		return lin_norm_stability(internals, comm_id);
	}

};


}
