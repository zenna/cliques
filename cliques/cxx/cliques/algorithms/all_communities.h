/* Copyright (c) Zenna Tavares-zennatavares@gmail.com, Michael Schuab 2010-2011 */
#pragma once

#include <vector>
#include <algorithm>


namespace cliques {

template<typename G>
void find_connected_communities(G &graph, std::vector<bool> allowed,
		std::vector<int> community,
		std::vector<std::vector<int> > & community_list) {
	bool any_allowed = false;
	bool breakout = false;
	int neigh_id = 0;
	int allowed_neigh = -1;

	if (community_list.size() % 10000 == 0) {
		cliques::output(community_list.size());
	}

	if (community.size() == 0) {
		for (unsigned int i = 0; i < allowed.size(); ++i) {
			if (allowed[i] == true) {
				//                std::cout << "start node" << i << std::endl;
				std::vector<int> start = { i };
				community_list.push_back(start);
				allowed[i] = false;
				find_connected_communities(graph, allowed, start,
						community_list);
			}
		}
	}

	for (unsigned int i = 0; i < community.size(); ++i) {
		typename G::Node n = graph.nodeFromId(community[i]);
		for (typename G::IncEdgeIt e(graph, n); e != lemon::INVALID; ++e) {
			typename G::Node neighbour = graph.oppositeNode(n, e);
			neigh_id = graph.id(neighbour);
			if (allowed[neigh_id] == true) {
				//                std::cout << "adding node" << neigh_id << std::endl;
				allowed_neigh = neigh_id;
				breakout = true;
				break;
			}
		}
		if (breakout) {
			break;
		}

	}
	if (allowed_neigh != -1) {
		community.push_back(allowed_neigh);
		allowed[allowed_neigh] = false;
		community_list.push_back(community);
		find_connected_communities(graph, allowed, community, community_list);

		community.pop_back();
		//        std::cout << "last removed" << last_removed << std::endl;
		find_connected_communities(graph, allowed, community, community_list);
	}

	for (unsigned int i = 0; i < allowed.size(); ++i) {
		if (allowed[i] == true) {
			//            std::cout << "allowed" << i << std::endl;
			any_allowed = true;
			break;
		}
	}

	if (any_allowed == false || (any_allowed == true && allowed_neigh == -1)) {
		return;
	}
}

/**
 @brief  Find all connected subgraphs of a graph

 This enumerates all the connected subgraphs of a graph
 */
template<typename G>
std::vector<std::vector<int> > find_connected_communities(G &graph) {
	std::vector<std::vector<int> > community_list;
	std::vector<bool> allowed(lemon::countNodes(graph), true);
	std::vector<int> community;
	find_connected_communities(graph, allowed, community, community_list);
	return community_list;
}


}
