/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_KERNIGHANLIN_H
#define CLIQUES_KERNIGHANLIN_H

#include <vector>
#include <map>
#include <limits>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#include <cliques/graphhelpers.h>
#include <cliques/algorithms/louvain.h>

namespace cliques
{
/**
 @brief  Refines a given partition using Kernighan Lin style movements

 This is a local algorithm which is based on modularity maximation outlined
 in Networks An Introduction, Newman page 373.  This in turn is based on
 an algorithm for finding graph bisection by Kernighan and Lin.

 It must take an initial partition, which may be random, but is typically the
 output of another optimisation such as louvain, upon which local optimisations
 are made.

 It works by looking at all pairs of nodes, and finding how the transition of node1
 into the community of n2 will change the quality (e.g. stability).  Over all
 pairs it executes the single best transition.  This process is repeated
 N times, while not moving nodes that have already been moved.
 Over the entire history of N transitions, it then finds the highest quality
 partition.  The algorithm is then recursed with the partition acting as the
 original partition as long as the quality continues to improve.

 Unlike louvain, this algorithm can make transitions which result in losses
 in quality in short run.

 @param[in]  graph     graph to find partition of
 @param[in]  compute_quality     partition quality function object
 @param[out]  std::vector<partitions>     optimal partitions, last in vector is overall best
 */
template<typename P, typename T, typename W, typename QF, typename QFDIFF>
double refine_partition_kernighan_lin(
		T &graph, W &weights, QF compute_quality,
		QFDIFF compute_quality_diff,
		P input_partition, P output_partition) {
	typedef typename T::Node Node;
	typedef typename T::NodeIt NodeIt;
	typedef lemon::RangeMap<double> range_map;

	Internals internals(graph, weights);
	int num_nodes = lemon::countNodes(graph);
	P partition = input_partition;
	P buffer_partition = partition;

	double minimum_improve = 0.000001;
	double best_quality = current_quality = compute_quality(internals);
	std::set<node> moved_nodes; //TODO change to unordered_set

	//TODO: Share internals across algorithms
	//TODO: Compute internals from a partition
	//TODO: Compute quality from a partition

	for (int i=0; i<num_nodes; ++i) {
		double best_gain = -std::numeric_limits<float>::max();
		Node node_to_move;
		unsigned int comm_to_move_to = n2_comm_id;

		for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
			// Don't move already moved nodes
			if (moved_nodes.find(n1) == moved_nodes.end()) {
				continue;
			}
			for (NodeIt n2(graph); n2 != lemon::INVALID; ++n2) {
				int n1_comm_id = partition.find_set(n1);
				int n2_comm_id = partition.find_set(n2);
				if (n1_comm_id == n2_comm_id) {
					continue;
				}
				isolate_and_update_internals(graph, weights, node, internals,
						partition, comm_id);
				double isolation_loss = compute_quality_diff(internals,
						n1_comm_id, node_id);
				double gain = compute_quality_diff(internals,
						n2_comm_id, node_id);
				double absolute_gain = gain - isolation_loss;
				// Put node back - we don't want to move yet
				insert_and_update_internals(graph, weights, n1, internals,
						partition, n1_comm_id);
				if (absolute_gain > best_gain) {
					best_gain = absolute_gain;
					node_to_move = n1;
					comm_to_move_to = n2_comm_id;
				}
			}
		}

		isolate_and_update_internals(graph, weights, node_to_move, internals,
				partition, node_to_move_comm);
		insert_and_update_internals(graph, weights, node_to_move, internals,
				partition, comm_to_move_to);
		moved_nodes.insert(node_to_move);

		current_quality = current_quality + absolute_gain;

		if (current_quality > best_quality) {
			buffer_partition = partition;
			best_quality = current_quality;
		}
	}

	if (best_quality > 0) {
		return refine_partition_kernighan_lin();
	}
	else {
		return partition;
	}
}

}

#endif
