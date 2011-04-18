/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_SPACE_H
#define CLIQUES_SPACE_H

#include <cliques/graphhelpers.h>

namespace cliques {

/**
 @brief  Finds neighbours of a partition

 */
template<typename P>
void find_neighbours(
		P partition,
		boost::unordered_set<P, cliques::partition_hash,
				cliques::partition_equal> &neighbour_partitions) {
	P temp_partition = partition;

	for (int node = 0; node < partition.element_count(); ++node) {
		int initial_set = temp_partition.find_set(node);
		temp_partition.isolate_node(node);

		for (typename P::PartIterator pitr = temp_partition.begin(); pitr
				!= partition.end(); ++pitr) {

			if (*pitr != initial_set) {
				temp_partition.add_node_to_set(node, *pitr);
				if (temp_partition != partition &&
						cliques::is_partition_connected(temp_partition)) {
					neighbour_partitions.insert(temp_partition);
				}
			}
		}
		temp_partition.add_node_to_set(node, initial_set);
	}
}

/**
 @brief  Convert a set of partitions into a graph

 Given a set of partitions (perhaps found with cliques::find_connected_partitions),
 this function creates a graph with nodes representing partitions and edges created
 according to a moveset rule.  The default rule is the single node moveset
 see cliques::find_neighbours

 @param[in] all_partitons reference to unordered set of partitions
 @param[out] space output graph representing the space
 */
template<typename G, typename P>
void create_space(
		boost::unordered_set<P, cliques::partition_hash,
				cliques::partition_equal> &all_partitions, G &space) {

	typedef typename G::Node Node;
	typedef typename G::Edge Edge;
	typedef boost::unordered_set<P, cliques::partition_hash,
			cliques::partition_equal> partition_set;
	typedef typename boost::unordered_set<P, cliques::partition_hash,
			cliques::partition_equal>::iterator partition_set_itr;

	boost::unordered_map<P, Node, cliques::partition_hash,
			cliques::partition_equal> partition_to_node;

	// Create node in space for each partition
	for (partition_set_itr itr = all_partitions.begin(); itr
			!= all_partitions.end(); ++itr) {
		Node temp_node = space.addNode();
		//partition_to_graph[temp_node] =
	}

	for (partition_set_itr itr = all_partitions.begin(); itr
			!= all_partitions.end(); ++itr) {

		Node current_node = partition_to_node[*itr];
		partition_set neighbour_partitions;
		cliques::find_neighbours(*itr, neighbour_partitions);

		for (partition_set_itr itr = neighbour_partitions.begin(); itr
				!= all_partitions.end(); ++itr) {
			Node neighbour_node = partition_to_node[*itr];
			// CHECK NOT SAME PARTITION
			Edge e = lemon::findEdge(space, current_node, neighbour_node);
			if (e == lemon::INVALID) {
				Edge e = space.addEdge(current_node, neighbour_node);
			}
		}
	}
}

template<typename G>
void create_disconnectivity_graph(G graph_landscape, G graph_dg) {
}

}

#endif
