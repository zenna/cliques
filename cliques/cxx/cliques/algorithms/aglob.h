/* Copyright (c) Zenna Tavares-zennatavares@gmail.com, Michael Schuab 2010-2011 */
#pragma once

#include <vector>
#include <limits>
#include <map>

namespace cliques {

template<typename G>
cliques::VectorPartition community_to_partition(G &graph,
		std::vector<int> comm, int outside_comm_set_id) {
	int num_nodes = lemon::countNodes(graph);
	cliques::VectorPartition partition(num_nodes, outside_comm_set_id);
	for (auto node = comm.begin(); node != comm.end(); ++node) {
		partition.add_node_to_set(*node, 1);
	}
	//    partition.normalise_ids();
	return partition;
}

/**
 @brief  Will removal of node break a community
 */
template<typename G, typename NO>
bool will_removal_break_community(G &graph, std::vector<int> comm,
		NO &node_to_remove) {
	cliques::VectorPartition partition = community_to_partition(graph, comm, 0);
	return cliques::will_move_break_partition(graph, partition, node_to_remove);
}

template<typename G>
std::vector<std::vector<int> > find_community_neighbours(G &graph,
		std::vector<int> comm) {

	std::set<int> forbidden;
	std::vector<std::vector<int> > neighbours;

	// Handle special case of single node community
	if (comm.size() == 1) {
		forbidden.insert(comm[0]);
	}
	for (unsigned int i = 0; i < comm.size(); ++i) { // Better if I could just iterate thro
		int node_id = comm[i];
		auto node = graph.nodeFromId(node_id);

		// If removal won't split community
		if (forbidden.find(node_id) == forbidden.end()
				&& will_removal_break_community(graph, comm, node) == false) {
			std::vector<int> contracted_comm = comm;
			auto location = std::find(contracted_comm.begin(),
					contracted_comm.end(), node_id);
			contracted_comm.erase(location);
			forbidden.insert(node_id);
			std::sort(contracted_comm.begin(), contracted_comm.end());
			//		    cliques::output("WHAAC");
			//		    cliques::print_collection(contracted_comm);
			neighbours.push_back(contracted_comm);
		}

		for (typename G::IncEdgeIt e(graph, node); e != lemon::INVALID; ++e) {
			auto neigh_node = graph.oppositeNode(node, e);
			int neigh_node_id = graph.id(neigh_node);

			// If node is on periphery
			if (std::find(comm.begin(), comm.end(), neigh_node_id)
					== comm.end() && forbidden.find(neigh_node_id)
					== forbidden.end()) {
				forbidden.insert(neigh_node_id);
				std::vector<int> new_comm = comm;
				new_comm.push_back(neigh_node_id);
				std::sort(new_comm.begin(), new_comm.end());

				neighbours.push_back(new_comm);//neighbours.push_back()
			}
		}
	}
	return neighbours;
}

/**
 @brief  Find maxima

 */
template<typename G, typename M, typename QF>
bool is_community_maxima(G &graph, M &weights, std::vector<int> comm,
		QF compute_quality, double &stability) {
	cliques::VectorPartition comm_partition = community_to_partition(graph,
			comm, 0);
	auto internals = cliques::gen(compute_quality, graph, weights,
			comm_partition);
	double current_quality = compute_quality(internals, 1, 2.3);
	double best_quality = -std::numeric_limits<float>::max();
	stability = current_quality;

	auto neighs = find_community_neighbours(graph, comm);
	for (auto itr = neighs.begin(); itr != neighs.end(); ++itr) {
		cliques::VectorPartition neigh_partition = community_to_partition(
				graph, *itr, 0);
		auto neigh_internals = cliques::gen(compute_quality, graph, weights,
				neigh_partition);
		double neigh_quality = compute_quality(neigh_internals, 1, 2.3);

		if (neigh_quality > best_quality) {
			best_quality = neigh_quality;
		}
	}

	// Neutrality?
	if (current_quality > best_quality) {
		return true;
	} else {
		return false;
	}

}

/**
 @brief  Optimise for maxima

 */
template<typename G, typename M, typename QF>
std::set<std::vector<int> > find_optimal_communities_huxley(G &graph,
		M &weights, QF &compute_quality, double time,
		std::vector<std::vector<int> > communities) {

	std::set<std::vector<int> > all_maxima;
	std::vector<int> comm = { 0 };
	int num_iterations = 1000;
	std::vector<int> buffer(lemon::countNodes(graph), 0);
	std::map<std::vector<int>, int> maxima_to_seen_count;
	for (int i = 0; i < num_iterations; ++i) {
		//		cliques::print_collection(comm);
		cliques::VectorPartition p = cliques::community_to_partition(graph,
				comm, 0);
		auto neighs = find_community_neighbours(graph, comm);
		double current_quality = compute_quality(p, 1, time);
		double real_quality = current_quality;
		for (auto maximum = maxima_to_seen_count.begin(); maximum
				!= maxima_to_seen_count.end(); ++maximum) {
			int dist_to_maximum = cliques::find_community_dist(graph, weights,
					comm, maximum->first, buffer);
			//			cliques::output("dist between following two is", dist_to_maximum);
			//			cliques::print_collection(comm);
			//			cliques::print_collection(maximum->first);
			current_quality -= discrete_gauss_kernel(dist_to_maximum, 0.5)
					* maximum->second;
		}

		std::vector<double> neigh_qualities;
		//		cliques::output("current quality:", current_quality);

		// Compute neighbour qualities then bias by history dependent maxima filling
		bool is_real_maximum = true;
		for (auto neigh = neighs.begin(); neigh != neighs.end(); ++neigh) {
			//			cliques::print_collection(*neigh);

			cliques::VectorPartition p = cliques::community_to_partition(graph,
					*neigh, 0);
			double neigh_quality = compute_quality(p, 1, time);

			if (neigh_quality > real_quality) {
				is_real_maximum = false;
			}
			for (auto maximum = maxima_to_seen_count.begin(); maximum
					!= maxima_to_seen_count.end(); ++maximum) {
				int dist_to_maximum = cliques::find_community_dist(graph,
						weights, *neigh, maximum->first, buffer);

				neigh_quality -= discrete_gauss_kernel(dist_to_maximum, 0.5)
						* maximum->second;
			}
			//			cliques::output("neigh quality:", neigh_quality);

			neigh_qualities.push_back(neigh_quality);
		}

		// Save if I am maxima on the unmodified landscape
		if (is_real_maximum == true) {
			all_maxima.insert(comm);
		}

		// Move to neighbour with probability dependent on difference in stability
		double best_quality_diff = -std::numeric_limits<double>::max();
		int j = 0, best_neighbour = -1;
		for (auto neigh_quality = neigh_qualities.begin(); neigh_quality
				!= neigh_qualities.end(); ++neigh_quality) {
			double quality_diff = *neigh_quality - current_quality;
			if (quality_diff > best_quality_diff) {
				best_quality_diff = quality_diff;
				best_neighbour = j;
			}
			++j;
		}
		// If not maxima in meta landscape
		if (best_quality_diff > 0.0) {
			comm = neighs[best_neighbour];
			cliques::output("moving", i);
			//			cliques::print_collection(comm);
		} else {
			cliques::output("meta maxima", i);
			maxima_to_seen_count[comm]++;
		}
	}

	return all_maxima;
}

template<typename G, typename M, typename QF>
std::set<std::vector<int> > find_optimal_communities_huxley(G &graph,
		M &weights, QF &compute_quality, double time) {
	std::vector<std::vector<int> > communities;
	return find_optimal_communities_huxley(graph, weights, compute_quality,
			time, communities);
}

}
