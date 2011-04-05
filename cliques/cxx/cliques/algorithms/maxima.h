/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010- 2011 */
#ifndef CLIQUES_MAXIMA_H
#define CLIQUES_MAXIMA_H

#include <cliques/structures/partition.h>
#include <limits>
#include <vector>
#include <lemon/concepts/graph.h>
#include <set>

namespace cliques
{
/**
@brief  Finds the maxima of a landscape graph through steepest ascent

This is an exhaustive search to find maximal points on a landscape
when a landscape is represented as a graph

@param[in] graph the input graph
@param[in] stabilities pointer (array) of node 'energies'
@tparam G graph type
*/
template <typename G>
std::set<int> find_maxima(G &graph, float *stabilities) {
    std::set<int> maxima;
    int num_iterations = 0;
	for (typename G::NodeIt n(graph); n != lemon::INVALID; ++n) {
		// Check part_itr hasn't already been visited
		if (num_iterations % 100 == 0) {
			//std::cout << "the number of iterations is" << num_iterations << std::endl;
			//std::cout << "so far number of maxima is " << maxima.size() << std::endl;
			//std::cout << "size of visited is " << visited.size() << std::endl;
			//std::cout << "num skipped " << num_skipped << std::endl;
		}
		typename G::Node best_neighbour = n;
		while (1) {
			bool has_improved = false;
			float best_score = stabilities[graph.id(best_neighbour)];
			for (typename G::OutArcIt a(graph, best_neighbour); a != lemon::INVALID; ++a) {
				if (best_score < stabilities[graph.id(graph.target(a))]) {
					best_score = stabilities[graph.id(graph.target(a))];
					best_neighbour = graph.target(a);
					has_improved = true;
				}
			}

			if (has_improved == false) {
				maxima.insert(graph.id(best_neighbour));
				break;
			}
		}
		num_iterations++;
	}

	return maxima;
};

}

#endif
