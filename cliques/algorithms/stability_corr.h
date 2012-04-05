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
#include <cliques/algorithms/internals/linearised_internals_corr.h>

namespace cliques {
/**
 @brief  Functor for finding normalised linearised correlation stability of weighted graph
 */
struct find_linearised_normalised_corr_stability {
	double markov_time;

	find_linearised_normalised_corr_stability(double markov_time) :
		markov_time(markov_time) {
	}

	template<typename G, typename P, typename W>
	double operator ()(G &graph, P &partition, W &weights,P &partition_init, std::vector<double> null_model) {
		cliques::LinearisedInternalsCorr internals(graph, weights, partition, partition_init,null_model);
		return (*this)(internals);
	}

	template<typename I>
	double operator ()(I &internals) {
		double q = 0;
		double q2 = 0;
		// first part
		int size = internals.null_model.size();
		for (int i = 0; i < size; i++) {
			// due to corr this first part is not equal 1-t any more
			q += 1 / (1 - internals.null_model[i]);
		}

		// second part
		size = internals.comm_loss.size();
		for (int i = 0; i < size; i++) {
			// main linear part, identical for cov. stability..
			if (internals.comm_loss[i] > 0) {
				q2 += markov_time * double(internals.comm_w_in[i])
						- (double(internals.comm_loss[i])
								* double(internals.comm_loss[i]));
			}
//			cliques::output("number",i,"in", internals.comm_w_in[i], "tot",internals.comm_loss[i], "nodew", internals.node_to_w[i],"q2", q2);
		}

		return q * (1.0 - markov_time) + q2;
	}

};

/**
 @brief  Functor for finding stability gain (normalised Laplacian) with for weighted graph
 */
struct linearised_normalised_corr_stability_gain {
	double markov_time;

	linearised_normalised_corr_stability_gain(double mtime) :
		markov_time(mtime) {
	}

	template<typename I>
	double operator ()(I &internals, int comm_id_neighbour, int node_id) {
		// loss factor so far...
		double comm_loss = internals.comm_loss[comm_id_neighbour];
		// gain from node
		double w_node_to_comm =
				internals.node_weight_to_communities[comm_id_neighbour];
		// weighted degree of node
		double w_deg_node = internals.node_to_w[node_id];
		return (markov_time * w_node_to_comm - comm_loss * w_deg_node) * 2;
	}
};

}
