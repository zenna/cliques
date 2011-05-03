/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010- 2011 */
#ifndef CLIQUES_HSG_H
#define CLIQUES_HSG_H

#include <vector>
#include <map>
#include <set>

#include <cliques/algorithms/louvain.h>

#include <cliques/algorithms/stability.h>
#include <cliques/structures/disjointset.h>

#include <cliques/helpers.h>

namespace cliques {
/**
@brief  Create a hierarchical layers using stability at different time-scales

@param[in] graph The input graph
@param[in] markov_times The Markov Times to optimise over
@param[out] layer_stabilities stability values for each layer

@return vector of partitions of increasing coarseness
*/
template <typename P, typename G>
std::vector<P> create_hsg_layers(G &graph, std::vector<float> &markov_times,
		std::vector<float> &layer_stabilities) {

	std::vector<P> partitions;
	for (std::vector<float>::iterator itr = markov_times.begin(); itr != markov_times.end(); ++itr) {

		std::vector<float> current_time;
		current_time.push_back(*itr);
		P best_partition = cliques::find_optimal_partition_louvain<P>(graph, cliques::find_linearised_stability(current_time));
		partitions.push_back(best_partition);
		layer_stabilities.push_back(*itr);
	}
	return partitions;
}

/**
@brief  Create a hierarchical stability graph

@param[in] graph The input graph
@param[in] markov_times Time range to create graph over
@param[out] hsg The output graph
@param[out] positions Node positions
*/
template <typename G>
void create_hsg (G &graph, std::vector<float> &markov_times, G &hsg,
		std::map<int,std::vector<float> > positions) {
    std::vector<float> layer_stabilities;

    typedef cliques::DisjointSetForest<int> DjForest;
	typedef typename G::Node Node;
	typedef typename G::Edge Edge;
	std::vector<DjForest> layers = create_hsg_layers<DjForest>(graph, markov_times, layer_stabilities);

	std::map<Node,Node> graph_to_hsg;

	for (std::vector<DjForest>::iterator partition = layers.begin(); partition != layers.end(); ++partition) {
		cliques::print_partition(*partition);
		for (DjForest::PartIterator pitr = partition->begin(); pitr != partition->end(); ++pitr) {
			Node new_node = hsg.addNode();
			std::cout << "creating node " << hsg.id(new_node) << std::endl;

    		for (DjForest::NodeIterator nitr = pitr.begin(); nitr != pitr.end(); ++nitr ) {
    			//First layer is special case
    			Node node_in_graph = graph.nodeFromId(*nitr);
    			if (partition == layers.begin()) {
    				graph_to_hsg[node_in_graph] = new_node;
    				continue;
    			}
    			Node node_in_hsg = graph_to_hsg[node_in_graph];

                Edge e = lemon::findEdge(hsg,node_in_hsg,new_node);
                if (e == lemon::INVALID) {
                    Edge e = hsg.addEdge(node_in_hsg,new_node);
                    //std::cout << "joining " << hsg.id(node_in_hsg) << " - "<< hsg.id(new_node) << std::endl;
                    //std::cout << "Edges:" << lemon::countEdges(hsg) << std::endl;
                }
    			graph_to_hsg[node_in_graph] = new_node;
    		}
    	}
	}
}

}

#endif
