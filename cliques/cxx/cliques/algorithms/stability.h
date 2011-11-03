/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/helpers.h>
#include<cliques/helpers/math.h>

#include <cliques/algorithms/internals/linearised_internals.h>
#include <cliques/algorithms/internals/linearised_internals_comb.h>
#include <cliques/graphhelpers.h>
#include <cliques/structures/disjointset.h>

namespace cliques {
/**
 @brief  Functor for finding normalised linearised stability of weighted graph
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
		//				cliques::output("time", markov_time, "size", size);

		for (int i = 0; i < size; i++) {
			if (internals.comm_w_tot[i] > 0) {
				q += markov_time * double(internals.comm_w_in[i]
						/ internals.two_m) - ((double(internals.comm_w_tot[i])
						/ internals.two_m) * (double(internals.comm_w_tot[i])
						/ internals.two_m));
			}
			//						cliques::output("in", internals.comm_w_in[i], "tot",internals.comm_w_tot[i], "nodew", internals.node_to_w[i], "2m", internals.two_m,"q", q);
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

		if (internals.comm_w_tot[i] == internals.two_m) {
			return q;
		} else
			return q / ((internals.comm_w_tot[i] / internals.two_m));// /(1-internals.comm_w_tot[i] / internals.two_m));
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

/**
 @brief  Functor for finding full normalised stability of weighted graph
 */
struct find_full_normalised_stability {

	lemon::SmartGraph exp_graph;
	lemon::SmartGraph::EdgeMap<double> exp_graph_weight;
	std::vector<double> minus_t_D_inv_L;
	std::vector<double> node_weighted_degree;
	find_linearised_normalised_stability lin_norm_stability;
	double m_time = -1;
	double threshold = 1e-9;

	template<typename G, typename M>
	find_full_normalised_stability(G &graph, M &weights, double thres) :
		lin_norm_stability(1) {

		// get threshold
		threshold = thres;

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

		// call expokit
		int N = node_weighted_degree.size();
		std::vector<double> exp_graph_vec = cliques::exp(minus_t_D_inv_L,
				markov_time, N);

		//		cliques::print_collection(exp_graph_vec);

		// create new graph out of matrix find_stability<lemon::SmartGraph>(partiton,time)
		lemon::SmartGraph::EdgeMap<double> exp_graph_weights(exp_graph);

		for (int i = 0; i < N; ++i) {
			for (int j = i; j < N; ++j) {
				double weight = exp_graph_vec[N * i + j]
						* node_weighted_degree[j];
				if (weight > threshold) {
					lemon::SmartGraph::Edge edge = exp_graph.addEdge(
							exp_graph.nodeFromId(i), exp_graph.nodeFromId(j));
					exp_graph_weights.set(edge, weight);
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

		// call expokit
		int N = node_weighted_degree.size();

		//        cliques::output("minus");
		//        cliques::print_collection(minus_t_D_inv_L, N);
		std::vector<double> exp_graph_vec = cliques::exp(minus_t_D_inv_L,
				markov_time, N);
		//        cliques::output("exp");
		//        cliques::print_collection(exp_graph_vec, N);

		// create new weight map out of matrix find_stability<lemon::SmartGraph>(partiton,time)


		if (markov_time != m_time) {
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
			exp_graph_weights(exp_graph);

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

//////////////////////////////////////////////////////////////////////////////
// Combinatorial Stability
//////////////////////////////////////////////////////////////////////////////
/**
 @brief  Functor for finding combinatorial linearised stability of weighted graph
 */
struct find_linearised_combinatorial_stability {
	double markov_time;

	find_linearised_combinatorial_stability(double markov_time) :
		markov_time(markov_time) {
	}

	template<typename G, typename P, typename W>
	double operator ()(G &graph, P &partition, W &weights, P &partition_init) {
		cliques::LinearisedInternalsComb internals(graph, weights, partition,
				partition_init);
		return (*this)(internals);
	}

	template<typename I>
	double operator ()(I &internals) {
		double q = 1.0 - markov_time;
		int size = internals.comm_tot_nodes.size();
		//				cliques::output("time", markov_time, "size", size);

		for (int i = 0; i < size; i++) {
			//			cliques::print_collection(internals.comm_tot_nodes);
			if (internals.comm_tot_nodes[i] > 0) {
				q += markov_time * double(internals.comm_w_in[i]
						/ internals.two_m)
						- ((double(internals.comm_tot_nodes[i])
								/ internals.num_nodes_init)
								* (double(internals.comm_tot_nodes[i])
										/ internals.num_nodes_init));
			}
			//						cliques::output("in", internals.comm_w_in[i], "tot",internals.comm_w_tot[i], "nodew", internals.node_to_w[i], "2m", internals.two_m,"q", q);
		}

		return q;
	}
};

/**
 @brief  Functor for finding stability gain (combinatorial Laplacian) with for weighted graph
 */
struct linearised_combinatorial_stability_gain {
	double markov_time;

	linearised_combinatorial_stability_gain(double mtime) :
		markov_time(mtime) {
	}

	double operator ()(double comm_tot_nodes, double w_node_to_comm,
			double two_m, double node_nr_nodes_init,
			unsigned int num_nodes_init) {
		return (markov_time * w_node_to_comm / two_m - (comm_tot_nodes
				/ num_nodes_init) * (node_nr_nodes_init / num_nodes_init));
	}

	template<typename I>
	double operator ()(I &internals, int comm_id_neighbour, int node_id) {
		double comm_tot_nodes = internals.comm_tot_nodes[comm_id_neighbour];
		double w_node_to_comm =
				internals.node_weight_to_communities[comm_id_neighbour];
		double node_nr_nodes_init = internals.node_to_nr_nodes_init[node_id];
		return (markov_time * w_node_to_comm / internals.two_m
				- (comm_tot_nodes / internals.num_nodes_init)
						* (node_nr_nodes_init / internals.num_nodes_init)) * 2;
	}
};

}
