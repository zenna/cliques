#ifndef CLIQUES_SPACE_H
#define CLIQUES_SPACE_H

namespace cliques {

/*template <typename G>
void create_landscape(umap &partition_map, G &graph_landscape) {
	for (umap::const_iterator part_itr = partition_map.begin(); part_itr != partition_map.end(); ++part_itr ) {

		if (num_iterations % 1000 == 0) {
			std::cout << "the number of iterations is" << num_iterations << std::endl;
		}
		uset partition_neighbours;
		find_neighbours(part_itr, partition_map, partition_neighbours);
		Node node_u = g_landscape.nodeFromId(part_itr->second->key);
		for (uset::iterator itr = partition_neighbours.begin(); itr != partition_neighbours.end();++itr) {
			Node node_v = g_landscape.nodeFromId((*itr)->second->key);
			Edge e = findEdge(g_landscape,node_u,node_v);
			if (e == INVALID) {
				Edge e = g_landscape.addEdge(node_u,node_v);
			}
		}
		num_iterations++;
	}
}*/

/**
@brief  Finds neighbours of a partition

*/
template <typename P>
void find_neighbours(P partition, boost::unordered_set<P, cliques::partition_hash, cliques::partition_equal> &neighbour_partitions) {
	int alpha;
	/*for each node
		Find its set -> set
		isolate node - change its parent to new parent / remove all children to first child
		for all sets other than set
			add node to set
			insert into neighbour partition
			undo
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
template <typename G, typename P>
void create_space(boost::unordered_set<P, cliques::partition_hash, cliques::partition_equal> &all_partitions, G &space) {

	typedef typename G::Node Node;
	typedef typename G::Edge Edge;
	boost::unordered_map<P,node,cliques::partition_hash, cliques::partition_equal> partition_to_node;

    for (typename boost::unordered_set<P,cliques::partition_hash, cliques::partition_equal>::iterator itr = all_partitions.begin();
     		itr != all_partitions.end(); ++itr) {
    	Node temp_node = space.addNode();
    	//partition_to_graph[temp_node] =
    }

    for (typename boost::unordered_set<P,cliques::partition_hash, cliques::partition_equal>::iterator itr = all_partitions.begin();
     		itr != all_partitions.end(); ++itr) {

    	Node current_node = partition_to_node[*itr];
    	boost::unordered_set<P,cliques::partition_hash, cliques::partition_equal> neighbour_partitions;
    	cliques::find_neighbours(*itr, neighbour_partitions);

    	for (typename boost::unordered_set<P,cliques::partition_hash, cliques::partition_equal>::iterator itr = neighbour_partitions.begin();
    	     		itr != all_partitions.end(); ++itr) {
    		Node neighbour_node = partition_to_node[*itr];
    		// CHECK NOT SAME PARTITION
            Edge e = lemon::findEdge(space, current_node, neighbour_node);
            if (e == lemon::INVALID) {
                Edge e = space.addEdge(current_node, neighbour_node);
            }
    	}
    }
}

template <typename G>
void create_disconnectivity_graph(G graph_landscape, G graph_dg) {
}

}

#endif //CLIQUES_COMPLEXITY_H
