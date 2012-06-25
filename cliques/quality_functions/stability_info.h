/* Copyright (c) Michael Schaub - michael.schaub09@imperial.ac.uk, 2010-2011 */
#pragma once
#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/helpers/helpers.h>
#include <cliques/helpers/graphhelpers.h>
#include <cliques/quality_functions/internals/linearised_internals_info.h>

namespace clq {
/**
 @brief  Functor for finding mutual information stability of weighted graph
 */
struct find_mutual_information_stability {

	template<typename G, typename P, typename W>
	double operator ()(G &graph, P &partition, W &weights) {
		clq::LinearisedInternalsInfo internals(graph, weights, partition);
		return (*this)(internals);
	}

	template<typename I>
	double operator ()(I &internals) {
		double q2 = 0;

		// second part
		int size = internals.comm_w_in.size();
		for (int i = 0; i < size; i++) {
			// main linear part, identical for cov. stability..
			if (internals.comm_w_in[i] > 0) {
				q2 += double(internals.comm_w_in[i]);
			}
			//			clq::output("number", i, "in", internals.comm_w_in[i], "q2", q2);
		}

		return q2;
	}

};

/**
 @brief  Functor for finding stability gain (normalised Laplacian) with for weighted graph
 */
struct mutual_information_stability_gain {

	template<typename I>
	double operator ()(I &internals, int comm_id_neighbour, int node_id) {
		// gain from node
		double w_node_to_comm =
				internals.node_weight_to_communities[comm_id_neighbour];
		return (w_node_to_comm) * 2;
	}
};


}
