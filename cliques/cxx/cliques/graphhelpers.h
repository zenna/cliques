#ifndef CLIQUES_GRAPHHELPERS_H
#define CLIQUES_GRAPHHELPERS_H

#include <iostream>
#include <fstream>
#include <map>

namespace cliques {

template <typename G>

float A(G &graph, int node1_id, int node2_id) {
    typename G::Node n1 = graph.nodeFromId(node1_id);
    typename G::Node n2 = graph.nodeFromId(node2_id);

    for(typename G::IncEdgeIt e(graph, n1); e!= lemon::INVALID; ++e) {
        if (graph.runningNode(e) == n2) {
            return 1.0;
        }
    }
    return 0.0;
};

template <typename G, typename M, typename NO>
float find_weighted_degree(G &graph, M &weights, NO node) {
	float degree = 0.0;
	for (typename G::IncEdgeIt e(graph,node); e != lemon::INVALID; ++e) {
		degree = degree + weights[e];
	}
	return degree;
};

template <typename G, typename M>
float find_total_weight(G &graph, M &weights) {
	float total_weight = 0.0;
	for (typename G::EdgeIt e(graph); e != lemon::INVALID; ++e) {
		total_weight = total_weight + weights[e];
	}
	return total_weight;
}


template <typename G>
void read_edgelist(G &graph, std::string filename) {

    std::ifstream maxima_file(filename.c_str());
    std::string line;
    std::string maxima;

    if ( !maxima_file.is_open() ) {
        std::cout << "couldn't open file" << std::endl;
        exit(1);
    }

    typedef typename G::Node Node;
    std::map<int,Node> id_to_node;
    while (std::getline(maxima_file,line)) {
        std::stringstream lineStream(line);
        std::set<int> current_maxima;

        std::getline(lineStream,maxima,' ');
        int node1_id = atoi(maxima.c_str());
        std::getline(lineStream,maxima,' ');
        int node2_id = atoi(maxima.c_str());


        typename std::map<int,Node>::iterator itr = id_to_node.find(node1_id);
        Node node1, node2;

        if (itr == id_to_node.end()) {
        	node1 = graph.addNode();
        	id_to_node[node1_id] = node1;
        }
        else {
        	node1 = itr->second;
        }
        itr = id_to_node.find(node2_id);
        if (itr == id_to_node.end()) {
        	node2 = graph.addNode();
        	id_to_node[node2_id] = node2;
        }
        else {
        	node2 = itr->second;
        }
        graph.addEdge(node1,node2);
    }
    maxima_file.close();
}

template <typename G, typename E>
void read_edgelist_weighted(std::string filename, G &graph, E &weights) {
	std::ifstream maxima_file(filename.c_str());
	std::string line;
	std::string maxima;

	if ( !maxima_file.is_open() ) {
		std::cout << "couldn't open file" << std::endl;
		exit(1);
	}

	typedef typename G::Node Node;
	std::map<int,Node> id_to_node;
	while (std::getline(maxima_file,line)) {
		std::stringstream lineStream(line);
		std::set<int> current_maxima;

		std::getline(lineStream,maxima,'\t');
		int node1_id = atoi(maxima.c_str());
		std::getline(lineStream,maxima,'\t');
		int node2_id = atoi(maxima.c_str());
		std::getline(lineStream,maxima,'\t');
		int weight = atoi(maxima.c_str());

		typename std::map<int,Node>::iterator itr = id_to_node.find(node1_id);
		Node node1, node2;

		if (itr == id_to_node.end()) {
			node1 = graph.addNode();
			id_to_node[node1_id] = node1;
		}
		else {
			node1 = itr->second;
		}
		itr = id_to_node.find(node2_id);
		if (itr == id_to_node.end()) {
			node2 = graph.addNode();
			id_to_node[node2_id] = node2;
		}
		else {
			node2 = itr->second;
		}
		typename G::Edge edge = graph.addEdge(node1,node2);
		weights.set(edge,weight);
	}
	maxima_file.close();
}


}

#endif
