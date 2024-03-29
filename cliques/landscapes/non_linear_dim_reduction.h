/* Copyright (c) Zenna Tavares, zennatavares@gmail.com, 2010-2011 */
#pragma once

#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>

#include "armadillo"
#include <lemon/dijkstra.h>
#include <lemon/concepts/graph.h>

#include <cliques/helpers/helpers.h>
#include <cliques/helpers/edit_distance.h>
#include <cliques/helpers/make_graphs.h>

//Issues
// Sampling depends on order of nodes considered
// Seems unsatisfactory, eigenvectors (but not eigenvalues) different
// Overlapping nodes makes visualisation difficult
// Sol: Z-order
// Visualisation slow

// TODO: Extend to larger graphs
// Triangulate
// Parallelise

// TODO: Error
// Sampled error
// Compute reconstruction error
//KNN error, compute extent to which neighbours in low_d space are also neighbours in high_d space
// For subset of nodes
// Find k neighbours in low_d and high d
// Compute score of overlap

namespace clq {

/**
 @brief  Find the pairwise geodesic (shortest path) distances between subset of nodes
 */
template<typename G, typename M>
arma::mat convert_graph_to_geodesic_dists(G &graph, std::vector<
		typename G::Node> nodes, M &weights) {
	lemon::concepts::ReadWriteMap<typename G::Node, int> dist_map;
	int N = nodes.size();
	arma::mat X(N, N);
	X.zeros();

	int num = 0;
	int total = N * (N - 1) / 2;
	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			auto d = lemon::Dijkstra<G, M>(graph, weights);
			auto n1 = nodes[i];
			auto n2 = nodes[j];
			d.run(n1, n2);
			float dist = d.dist(n2);
			X(i, j) = dist;
			X(j, i) = dist;
			if (num % 1000000 == 0) {
				clq::output(graph.id(n1), graph.id(n2), dist, num, ":",
						total);
			}
			num++;
		}
	}
	return X;
}

/**
 @brief  Find the euclidean distance between two row vectors
 */
inline
double euclid_dist(arma::rowvec A, arma::rowvec B) {
	int length = A.n_elem;
	double result = 0.0;
	for (int i = 0; i < length; ++i) {
		double temp = A(i) - B(i);
		result += temp * temp;
	}
	return std::sqrt(result);
}

/**
 @brief  Find the pairwise distances between all row pairs of a matrix
 */
arma::mat euclid_pairwise_dists(arma::mat A) {
	int N = A.n_rows;
	arma::mat D(N, N);
	D.zeros();
	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			D(i, j) = euclid_dist(A.row(i), A.row(j));
			D(j, i) = D(i, j);
		}
	}
	return D;
}

/**
 @brief  Find the standard linear correlation coefficient
 */
double correlation_coeff(arma::mat &A, arma::mat &B) {
	double mean_A = arma::mean(arma::mean(A));
	double mean_B = arma::mean(arma::mean(B));

	double numerator_sum = 0.0;
	double denom_a_sum = 0.0;
	double denom_b_sum = 0.0;

	int n_rows = A.n_rows;
	int n_cols = A.n_cols;
	for (int i = 0; i < n_rows; ++i) {
		for (int j = 0; j < n_cols; ++j) {
			double one = A(i, j) - mean_A;
			double two = B(i, j) - mean_B;
			numerator_sum += one * two;
			denom_a_sum += one * one;
			denom_b_sum += two * two;
		}
	}
	return numerator_sum / std::sqrt(denom_a_sum * denom_b_sum);
}

/**
 @brief  Find the residual variance
 1 - R^2(D_m, D_y)
 D_m = ;
 D_y = matrix of euclidean distance in low_dim;
 */
double residual_variance(arma::mat &D_m, arma::mat &D_y) {
	return 1 - correlation_coeff(D_m, D_y);
}

/**
 @brief  Choose a random subset of nodes from a graph
 */
template<typename G>
std::vector<typename G::Node> randomly_choose_nodes(G &graph,
		int num_random_nodes) {
	int num_nodes = lemon::countNodes(graph);
	srand(time(NULL));

	std::vector<typename G::Node> landmark_nodes;
	for (int i = 0; i < num_random_nodes; ++i) {
		int rand_node_id = rand() % num_nodes;
		landmark_nodes.push_back(graph.nodeFromId(rand_node_id));
	}
	return landmark_nodes;
}

//[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
//view.js:1080 0 44 44
//view.js:114[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0]

template<typename G, typename M, typename C>
double find_community_dist(G &graph, M &weights, C commi, C commj, C &diff_itr) {
	int size1 = commi.size();
	int size2 = commj.size();

	std::sort(commi.begin(), commi.end());
	std::sort(commj.begin(), commj.end());

	auto end = std::set_intersection(commi.begin(), commi.end(),
			commj.begin(), commj.end(), diff_itr.begin());
	int intersection_size = int(end - diff_itr.begin());

	int difference_size = size1 + size2 - intersection_size * 2;
	float shortest_inter_comm_distance = 0.0;

	if (intersection_size == 0) {
		difference_size = size1 + size2 - 1;
		shortest_inter_comm_distance = std::numeric_limits<int>::max();
		for (auto n1 = commi.begin(); n1 != commi.end(); ++n1) {
			for (auto n2 = commj.begin(); n2 != commj.end(); ++n2) {
				auto d = lemon::Dijkstra<G, M>(graph, weights);
				typename G::Node node1 = graph.nodeFromId(*n1);
				typename G::Node node2 = graph.nodeFromId(*n2);
				d.run(node1, node2);
				double dist = d.dist(node2);
				if (shortest_inter_comm_distance > dist) {
					shortest_inter_comm_distance = dist;
				}
			}
		}
		difference_size = size1 + size2 + shortest_inter_comm_distance - 2;
	}

	float edit_distance = shortest_inter_comm_distance + difference_size;
	return edit_distance;
}

template<typename G, typename S>
arma::mat find_community_edit_dists(G &graph, S &communities) {
	//	typedef typename G::EdgeMap<float> EdgeMap;
	//	EdgeMap eights(graph);
	lemon::SmartGraph::EdgeMap<float> weights(graph);
	clq::make_weights_from_edges(graph, weights);

	int N = communities.size();
	arma::mat X(N, N);
	X.zeros();
	lemon::concepts::ReadWriteMap<typename G::Node, int> dist_map;
	int total = (N * (N - 1)) / 2;
	std::vector<int> diff_itr(lemon::countNodes(graph));

	int i = 0;
	int j = 0;
	int num = 0;
	for (auto commi = communities.begin(); commi != communities.end(); ++commi) {
		j = i + 1;
		for (auto commj = commi; ++commj != communities.end();) {
			float edit_distance = find_community_dist(graph, weights, *commi,
					*commj, diff_itr);
			X(i, j) = edit_distance;
			X(j, i) = edit_distance;
			if (num % 1000000 == 0) {
				clq::output(i, j, edit_distance, num, ":", total);
			}
			num++;
			j++;
		}
		i++;
	}
	return X;
}

/**
 @brief  Find the pairwise partition edit distances (equivalent to above for partitions)
 */
template<typename G, typename BM>
arma::mat find_edit_dists(G &graph, std::vector<typename G::Node> nodes,
		BM &map) {
	int N = nodes.size();
	arma::mat X(N, N);
	X.zeros();

	int num = 0;
	int total = N * (N - 1) / 2;
	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			auto n1 = nodes[i];
			auto n2 = nodes[j];
			auto p1 = map.right.at(n1);
			auto p2 = map.right.at(n2);
			clq::Hungarian hungarian(p1, p2);
			float edit_distance = float(hungarian.edit_distance());
			X(i, j) = edit_distance;
			X(j, i) = edit_distance;
			if (num % 1000000 == 0) {
				clq::output(graph.id(n1), graph.id(n2), edit_distance, num,
						":", total);
			}
			num++;
		}
	}
	return X;
}

/**
 @brief  Find the pairwise partition edit distances (equivalent to above for partitions)
 */
template<typename S>
arma::mat find_edit_dists(S &partitions) {
	int N = partitions.size();
	arma::mat X(N, N);
	X.zeros();

	int total = (N * (N - 1)) / 2;

	int i = 0;
	int j = 0;
	int num = 0;
	for (auto itr1 = partitions.begin(); itr1 != partitions.end(); ++itr1) {
		j = i + 1;
		for (auto itr2 = itr1; ++itr2 != partitions.end();) {
			auto p1 = *itr1;
			auto p2 = *itr2;
			clq::Hungarian hungarian(p1, p2);
			float edit_distance = float(hungarian.edit_distance());
			X(i, j) = edit_distance;
			X(j, i) = edit_distance;
			if (num % 1000000 == 0) {
				clq::output(i, j, edit_distance, num, ":", total);
			}
			num++;
			j++;
		}
		i++;
	}
	return X;
}

/**
 @brief  Find the pairwise partition edit distances (equivalent to above for partitions)
 */
template<typename S>
arma::mat find_edit_landmark_dists(S &partitions, S &landmarks) {
	int N = partitions.size();
	int N_l = landmarks.size();
	arma::mat D_l(N_l, N);
	D_l.zeros();

	int i = 0;
	for (auto l_itr = landmarks.begin(); l_itr != landmarks.end(); ++l_itr) {
		int j = 0;
		for (auto n_itr = landmarks.begin(); n_itr != landmarks.end(); ++n_itr) {
			auto p1 = *n_itr;
			auto p2 = *l_itr;
			clq::Hungarian hungarian(p1, p2);
			float edit_distance = float(hungarian.edit_distance());
			D_l(i, j) = edit_distance;
			++j;
		}
		++i;
	}
	return D_l;
}

///**
// @brief  Find the pairwise partition hasse diagram distances
// */
template<typename S>
arma::mat find_split_merge_dists(S &partitions) {
	// Get the number of partition (allowed)
	int N = partitions.size();
	// allocate pairwise distance matrix
	arma::mat X(N, N);
	X.zeros();
	// total number of nodes; needed for print statement
	int total = (N * (N - 1)) / 2;

	// pre-allocate some counters
	int i, j, num;
	i = j = num = 0;

	// loop over all pairs of partitions
	for (auto itr1 = partitions.begin(); itr1 != partitions.end(); ++itr1) {
		j = i + 1;
		for (auto itr2 = itr1; ++itr2 != partitions.end();) {

			// get the partitions for comparisons
			auto p1 = *itr1;
			auto p2 = *itr2;

			// get the difference of the partitions in terms of the number of communities
			int nr_comm_diff = p1.set_count() - p2.set_count();
			//clq::output("difference:", nr_comm_diff);

			// only partitions with difference one can be direct neighbours
			if (std::abs(nr_comm_diff) == 1
					&& are_partitions_connected_with_one_merge(p1, p2)) {
				X(j, i) = 1;
				X(i, j) = 1;
				num++;
//				clq::output("direct neighbour", num);
				if (num % 1000000 == 0) {
					clq::output(i, j, X(i, j), num, ":", total);
				}
			}


			j++;
		}
		i++;
	}

	arma::mat A(N, N);
	A= X;
	int pathlength = 1;
	// as long as we have not found all shortest paths continue.
	while (num != total) {
		A = A * X;
		pathlength++;
		for (int i = 0; i < N; ++i) {
			for (int j = i+1; j < N; ++j) {
				if (A(i, j) != 0 && X(i,j)==0) {
					X(i, j) = pathlength;
					X(j, i) = pathlength;
					num++;
				}
			}
		}
	}

	return X;
}

bool are_partitions_connected_with_one_merge(VectorPartition p1,
		VectorPartition p2) {
	// get difference in number of communities, only partitions with difference one can be direct neighbours
	int num_commi = p1.set_count();
	int num_commj = p2.set_count();
	if (std::abs(num_commi - num_commj) != 1) {
		return false;
	}
//	clq::output("comparison");
//	clq::print_partition_line(p1);
//	clq::print_partition_line(p2);
	VectorPartition greater_partition(p1.return_partition_vector());
	VectorPartition smaller_partition(p2.return_partition_vector());
	if (num_commi - num_commj == -1) {
		greater_partition = VectorPartition(p2.return_partition_vector());
		smaller_partition = VectorPartition(p1.return_partition_vector());
	}

	// initialise
	int difference_size = 0;

	// normalise ids
	greater_partition.normalise_ids();
	smaller_partition.normalise_ids();
	int maxi = greater_partition.set_count();
	int i = 0;
	while(i<maxi){
		// get the two communities
		std::vector<int> commi = greater_partition.get_nodes_from_set(i);
		std::vector<int> commj = smaller_partition.get_nodes_from_set(i);

//		clq::output("communities");
//		clq::print_collection(commi);
//		clq::print_collection(commj);

		// in here the difference is stored
		std::vector<int> difference(commi.size() + commj.size());
		std::vector<int>::iterator it;
		it = std::set_symmetric_difference(commi.begin(), commi.end(), commj.begin(),
				commj.end(), difference.begin());

		// check if communities match
		difference_size = int(it - difference.begin());

		// if they do not check difference
		if (difference_size != 0) {
			// if this is the only difference we got one more chance

				// cut difference to correct size
				difference = std::vector<int>(difference.begin(), it);
//				clq::print_collection(difference);

				// difference should be smaller then each community
				if (difference.size() >= commi.size() && difference.size() >= commj.size()){
					return false;
				}
				// check if commi was the bigger one, merge smaller communities
				if (commi.size() > commj.size()) {
					// merge communities in other partition
					for (int k = 0; k < int(commj.size()); ++k) {
						smaller_partition.add_node_to_set(commj[k], smaller_partition.find_set(difference[0]));
					}
					smaller_partition.normalise_ids();
				} else {
					// do it the other way around
					for (int k = 0; k < int(commi.size()); ++k) {
						greater_partition.add_node_to_set(commi[k], greater_partition.find_set(difference[0]));
					}
					greater_partition.normalise_ids();
				}

				for(int k = 0; k<greater_partition.element_count();++k){
					if(greater_partition.find_set(k)!=smaller_partition.find_set(k)){
						return false;
					}
				}


		}
	i++;
	}
	return true;
}

/**
 @brief  Embed pairwise distance matrix into another (lower) dimension
 //TODO deal with negative eigen values
 */
arma::mat embed_mds(arma::mat &X, int num_dimen) {
	int N = X.n_cols;
	arma::vec means = arma::zeros<arma::vec>(N);
	for (int i = 0; i < N; ++i) {
		means(i) = arma::mean(X.col(i));
	}
	double mean = arma::mean(means);
	arma::mat B_n(N, N);
	B_n.zeros();

	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			B_n(i, j) = -(X(i, j) - means(i) - means(j) + mean) / 2.0;
			B_n(j, i) = B_n(i, j);
		}
	}

	arma::vec eigvals;
	arma::mat eigvecs;
	eig_sym(eigvals, eigvecs, B_n);

	arma::uvec indices = arma::sort_index(eigvals, 1);
	arma::mat L(num_dimen, N);
	for (int i = 0; i < num_dimen; ++i) {
		for (int j = 0; j < N; ++j) {
			L(i, j) = eigvecs(j, indices(i)) * std::sqrt(eigvals(indices(i)));
		}
	}
	return L;
}

/**
 @brief  Embed pairwise distance matrix into another (lower) dimension
 */
arma::mat embed_landmark_mds(arma::mat &D_n, arma::mat &D_l, int num_dimen) {
	int N = D_n.n_cols;
	arma::vec means = arma::zeros<arma::vec>(N);
	for (int i = 0; i < N; ++i) {
		means(i) = arma::mean(D_n.col(i));
	}
	double mean = arma::mean(means);
	arma::mat B_n(N, N);
	B_n.zeros();

	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			B_n(i, j) = -(D_n(i, j) - means(i) - means(j) + mean) / 2.0;
			B_n(j, i) = B_n(i, j);
		}
	}

	arma::vec eigvals;
	arma::mat eigvecs;
	eig_sym(eigvals, eigvecs, B_n);

	arma::uvec indices = arma::sort_index(eigvals, 1);
	arma::mat L(num_dimen, N);
	for (int i = 0; i < num_dimen; ++i) {
		for (int j = 0; j < N; ++j) {
			L(i, j) = eigvecs(j, indices(i)) * std::sqrt(eigvals(indices(i)));
		}
	}

	arma::mat L_sharp(num_dimen, N);
	for (int i = 0; i < num_dimen; ++i) {
		for (int j = 0; j < N; ++j) {
			L_sharp(i, j) = eigvecs(j, indices(i)) / std::sqrt(eigvals(indices(
					i)));
		}
	}

	int N_all = D_l.n_cols;
	arma::mat X(num_dimen, N_all);
	for (int i = 0; i < N_all; ++i) {
		arma::vec second_term = D_l.col(i) - means;
		X.col(i) = -L_sharp * second_term / 2.0;
	}

	return X;
}

/**
 @brief  Embeds a graph in a lower dim space (convenience function)
 */
template<typename G, typename M>
arma::mat embed_graph(G& graph, M& weights, int num_dim) {
	typedef typename G::Node Node;
	typedef typename G::NodeIt NodeIt;

	// Use all nodes
	std::vector<Node> landmark_nodes;
	for (NodeIt n(graph); n != lemon::INVALID; ++n) {
		Node node = n;
		landmark_nodes.push_back(node);
	}

	clq::output("choosing nodes");
	auto X = clq::convert_graph_to_geodesic_dists(graph,
			landmark_nodes, weights);
	X.print("X");

	clq::output("finding embedding");
	auto L = clq::embed_mds(X, num_dim);
	L.print("L");

	clq::output("saving");
	arma::mat L_t = arma::trans(L);

	auto D_y = clq::euclid_pairwise_dists(L_t);
	clq::output("residual variance", clq::residual_variance(X, D_y));

	return L_t;
}

template<typename G, typename M, typename T>
void convert_dists_to_graph(G &graph, M &weights, arma::mat &dists, T threshold) {
	typedef typename G::Node Node;
	typedef typename G::Edge Edge;
	int N = dists.n_cols;
	for (int i = 0; i < N; ++i) {
		graph.addNode();
	}

	for (int i = 0; i < N - 1; ++i) {
		for (int j = i + 1; j < N; ++j) {
			if (dists(i, j) == threshold) {
				Node u = graph.nodeFromId(i);
				Node v = graph.nodeFromId(j);
				Edge new_edge = graph.addEdge(u, v);
				weights.set(new_edge, dists(i, j));
			}
		}
	}
}
}
