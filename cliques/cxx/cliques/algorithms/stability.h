/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

//#include <cliques/helpers.h>
#include<cliques/helpers/math.h>

#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/graphhelpers.h>
#include <cliques/structures/disjointset.h>

namespace cliques {
/**
 @brief  Functor for finding linearised stability of weighted graph
 */
struct find_linearised_normalised_stability {
	double markov_time;

	find_linearised_normalised_stability(double markov_time) :
		markov_time(markov_time) {
	}

	template<typename G, typename P, typename W>
	double operator ()(G &graph, P &partition, W &weights) {
		cliques::LinearisedInternals internals(graph, weights, partition);
		return (*this)(internals);
	}

	template<typename I>
	double operator ()(I &internals) {
		double q = 1.0 - markov_time;
		int size = internals.comm_w_tot.size();
		//		cliques::output("time", markov_time, "size", size);

		for (int i = 0; i < size; i++) {
			if (internals.comm_w_tot[i] > 0) {
				q += markov_time * double(internals.comm_w_in[i]
						/ internals.two_m) - ((double(internals.comm_w_tot[i])
						/ internals.two_m) * (double(internals.comm_w_tot[i])
						/ internals.two_m));
			}
			//			cliques::output("q", q);
		}

		return q;
	}

	template<typename I>
	double operator ()(I &internals, int comm_id) {
		double q = -1;
		int i = comm_id;
		if (internals.comm_w_tot[i] > 0) {
			//            cliques::output("internals", internals.comm_w_tot[i], internals.two_m, markov_time);
			q = internals.comm_w_tot[i] / internals.two_m * (1.0 - markov_time);
			q += markov_time * double(internals.comm_w_in[i] / internals.two_m)
					- ((double(internals.comm_w_tot[i]) / internals.two_m)
							* (double(internals.comm_w_tot[i])
									/ internals.two_m));
		} else {
			std::cout << "This community does not exist!!!" << std::endl;
		}
		//        cliques::output("Q", q);
		return q / internals.comm_w_tot[i] * internals.two_m;
	}

};

/**
 @brief  Functor for finding stability gain (normalised Laplacian) with for weighted graph
 */
struct linearised_normalised_stability_gain {
	double markov_time;

	linearised_normalised_stability_gain(double mtime) :
		markov_time(mtime) {
	}

	double operator ()(double tot_w_comm, double w_node_to_comm, double two_m,
			double w_deg_node) {
		return (markov_time * w_node_to_comm - tot_w_comm * w_deg_node / two_m);
	}

	template<typename I>
	double operator ()(I &internals, int comm_id_neighbour, int node_id) {
		double tot_w_comm = internals.comm_w_tot[comm_id_neighbour];
		double w_node_to_comm =
				internals.node_weight_to_communities[comm_id_neighbour];
		double w_deg_node = internals.node_to_w[node_id];
		return (markov_time * w_node_to_comm - tot_w_comm * w_deg_node
				/ internals.two_m) * 2 / internals.two_m;
	}
};

struct find_full_normalised_stability {
	double markov_time;
	find_linearised_normalised_stability lin_norm_stability;

	find_full_normalised_stability(double markov_time) :
		markov_time(markov_time), lin_norm_stability(1) {
	}

	template<typename G, typename P, typename W>
	double operator ()(G &graph, P &partition, W &weights) {

		typedef typename G::Node Node;
		typedef typename G::NodeIt NodeIt;
		typedef typename G::EdgeIt EdgeIt;
		/////////////////
		// read out graph into matrix
		//////////////////

		// number of nodes N
		int N = lemon::countNodes(graph);
		// -D^-1*L*t == t(B-I)
		std::vector<double> minus_t_D_inv_L(N * N, 0);
		std::vector<double> node_weighted_degree(N, 0);

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
			// t*B_uv
			minus_t_D_inv_L[node_id_v + N * node_id_u] += weight_uv
					/ node_weighted_degree[node_id_v];
			// t*B_vu
			minus_t_D_inv_L[node_id_u + N * node_id_v] += weight_uv
					/ node_weighted_degree[node_id_u];
		}

		// call expokit
		std::vector<double> exp_graph_vec = cliques::exp(minus_t_D_inv_L, markov_time,
				N);

//		cliques::print_collection(exp_graph_vec);

		// create new graph out of matrix
		G exp_graph;
		typename G::template EdgeMap<double>  exp_graph_weights(exp_graph);

		// reserve memory space for number of nodes
		exp_graph.reserveNode(N);
		// add nodes
		for (int i = 0; i < N; ++i) {
			exp_graph.addNode();
		}

		for (int i = 0; i < N; ++i) {
			for (int j = i + 1; j < N; ++j) {
				double weight = exp_graph_vec[N * i + j]
						* node_weighted_degree[j];
				if (weight > 0) {
					typename G::Edge edge = exp_graph.addEdge(exp_graph.nodeFromId(i), exp_graph.nodeFromId(j));
					exp_graph_weights.set(edge, weight);
				}
			}
		}

		cliques::LinearisedInternals internals(exp_graph, exp_graph_weights,
				partition);
		return lin_norm_stability(internals);
	}

};

}
