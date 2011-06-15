/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/graphhelpers.h>
#include <cliques/structures/disjointset.h>

namespace cliques {
/**
 @brief  Functor for finding stability of weighted graph
 */
struct find_weighted_linearised_stability {
	std::vector<double> &markov_times_;

	find_weighted_linearised_stability(std::vector<double> &markov_times) :
		markov_times_(markov_times) {
	}

	template<typename G, typename T2, typename W>
	void operator ()(G &graph, T2 &my_partition, W &weights,
			std::vector<double> &stabilities) {
		float two_m = 2.0 * cliques::find_total_weight(graph, weights);
		float first_term = 0.0;
		float second_term = 0.0;

		for (typename T2::PartIterator pitr = my_partition.begin(); pitr
				!= my_partition.end(); ++pitr) {
			for (typename T2::NodeIterator n1itr = pitr.begin(); n1itr
					!= pitr.end(); ++n1itr) {
				for (typename T2::NodeIterator n2itr = pitr.begin(); n2itr
						!= pitr.end(); ++n2itr) {
					//std::cout << "n1: " << *n1itr << " n2: " << *n2itr << std::endl;
					float k1 = cliques::find_weighted_degree(graph, weights,
							graph.nodeFromId(*n1itr));
					float k2 = cliques::find_weighted_degree(graph, weights,
							graph.nodeFromId(*n2itr));
					typename G::Edge edge = lemon::findEdge(graph,
							graph.nodeFromId(*n1itr), graph.nodeFromId(*n2itr));
					float A = weights[edge];
					first_term = first_term + k1 * k2;
					second_term = second_term + A;
				}
			}
		}
		first_term = first_term / (two_m * two_m);
		second_term = second_term / two_m;

		float R;
		for (std::vector<double>::iterator t = markov_times_.begin(); t
				!= markov_times_.end(); ++t) {
			R = 0.0;
			R = (1.0 - *t) - first_term + *t * second_term;
			stabilities.push_back(R);
		}
	}

	template<typename I>
	double operator ()(I &internals, int i = 0) {
		double markov_time = markov_times_[i];
		double q = 1.0 - markov_time;
		int size = internals.comm_w_tot.size();

		for (int i = 0; i < size; i++) {
			if (internals.comm_w_tot[i] > 0) {
				q += double(internals.comm_w_in[i]
						/ internals.two_m) - ((double(internals.comm_w_tot[i])
						/ internals.two_m) * (double(internals.comm_w_tot[i])
						/ internals.two_m));
			}
			q = q * markov_time;
		}

		return q;
	}
};

/**
 @brief  Functor for finding stability gain for weighted graph
 */
struct linearised_stability_gain_louvain {
	double markov_time;

	linearised_stability_gain_louvain(double mtime) :
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
		return (markov_time *w_node_to_comm - tot_w_comm * w_deg_node
				/ internals.two_m)*2/internals.two_m;
	}
};

}
