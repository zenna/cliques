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

#include <cliques/helpers.h>
#include <cliques/helpers/edit_distance.h>
#include <cliques/structures/make_graphs.h>

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

namespace cliques {

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

/**
 @brief  Find the pairwise geodesic (shortest path) distances between subset of nodes
 */
template<typename G, typename M>
arma::mat find_geodesic_dists(G &graph, std::vector<typename G::Node> nodes,
        M &weights) {
    lemon::concepts::ReadWriteMap<typename G::Node, int> dist_map;
    int N = nodes.size();
    arma::mat X(N, N);
    X.zeros();

    int num = 0;
    int total = N * (N - 1) / 2;
    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            auto d = lemon::Dijkstra<G,M>(graph, weights);
            auto n1 = nodes[i];
            auto n2 = nodes[j];
            d.run(n1, n2);
            float dist = d.dist(n2);
            X(i,j) = dist;
            X(j,i) = dist;
            if (num  % 1000000 == 0) {
                cliques::output(graph.id(n1), graph.id(n2), dist, num, ":",total);
            }
            num++;
        }
    }
    return X;
}

//[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
//view.js:1080 0 44 44
//view.js:114[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0]

template <typename G, typename M, typename C>
float find_community_dist(G &graph, M &weights, C comm1, C comm2, C &diff_itr) {
    int size1 = comm1.size();
    int size2 = comm2.size();

    std::sort(comm1.begin(),comm1.end());
    std::sort(comm2.begin(),comm2.end());

	auto end = std::set_symmetric_difference(comm1.begin(),
			comm1.end(), comm2.begin(),
			comm2.end(), diff_itr.begin());
	int difference_size = int(end - diff_itr.begin());

    float shortest_inter_comm_distance = 0.0;

    if (difference_size == 0) {
    	difference_size = size1 + size2 - 1;
    	shortest_inter_comm_distance = std::numeric_limits<int>::max();
        for (auto n1 = comm1.begin(); n1 != comm1.end(); ++n1) {
        	for (auto n2 = comm2.begin(); n2 != comm2.end(); ++n2) {
                auto d = lemon::Dijkstra<lemon::SmartGraph,lemon::SmartGraph::EdgeMap<float> >(graph, weights);
                typename G::Node node1 = graph.nodeFromId(*n1);
                typename G::Node node2 = graph.nodeFromId(*n2);
                d.run(node1, node2);
                float dist = d.dist(node2);
                if (shortest_inter_comm_distance > dist) {
                	shortest_inter_comm_distance = dist;
                }
        	}
        }
    }

    float edit_distance = shortest_inter_comm_distance + difference_size;
    return edit_distance;
}

template <typename G, typename S>
arma::mat find_community_edit_dists(G &graph, S &communities) {
//	typedef typename G::EdgeMap<float> EdgeMap;
//	EdgeMap eights(graph);
	lemon::SmartGraph::EdgeMap<float> weights(graph);
	cliques::make_weights_from_edges(graph, weights);

	int N = communities.size();
	arma::mat X(N,N);
	X.zeros();
    lemon::concepts::ReadWriteMap<typename G::Node, int> dist_map;
	int total = (N * (N - 1)) / 2;
    std::vector<int> diff_itr(lemon::countNodes(graph));

    int i=0;
    int j=0;
    int num = 0;
    for (auto comm1 = communities.begin(); comm1 != communities.end(); ++comm1) {
    	j = i+1;
        for (auto comm2 = comm1; ++comm2 != communities.end();) {
//            int size1 = comm1->size();
//            int size2 = comm2->size();
//
//            std::sort(comm1->begin(),comm1->end());
//            std::sort(comm2->begin(),comm2->end());
//
//			auto end = std::set_difference(comm1->begin(),
//					comm1->end(), comm2->begin(),
//					comm2->end(), diff_itr.begin());
//			int difference_size = int(end - diff_itr.begin());
//
//            float shortest_inter_comm_distance = 0.0;
//
//            if (difference_size == 0) {
//            	difference_size = size1 + size2 - 1;
//            	shortest_inter_comm_distance = std::numeric_limits<int>::max();
//                for (auto n1 = comm1->begin(); n1 != comm1->end(); ++n1) {
//                	for (auto n2 = comm2->begin(); n2 != comm2->end(); ++n2) {
//                        auto d = lemon::Dijkstra<lemon::SmartGraph,lemon::SmartGraph::EdgeMap<float> >(graph, weights);
//                        typename G::Node node1 = graph.nodeFromId(*n1);
//                        typename G::Node node2 = graph.nodeFromId(*n2);
//                        d.run(node1, node2);
//                        float dist = d.dist(node2);
//                        if (shortest_inter_comm_distance > dist) {
//                        	shortest_inter_comm_distance = dist;
//                        }
//                	}
//                }
//            }
//
//            float edit_distance = shortest_inter_comm_distance + difference_size;
            float edit_distance = find_community_dist(graph,weights,*comm1, *comm2, diff_itr);
            X(i,j) = edit_distance;
            X(j,i) = edit_distance;
            if (num  % 1000000 == 0) {
                cliques::output(i,j, edit_distance, num, ":", total);
            }
            num++;
            j++;
        }i++;
    }
    return X;
}

/**
 @brief  Find the pairwise partition edit distances (equivalent to above for partitions)
 */
template<typename G, typename BM>
arma::mat find_edit_dists(G &graph, std::vector<typename G::Node> nodes, BM &map) {
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
            cliques::Hungarian hungarian(p1,p2);
            float edit_distance = float(hungarian.edit_distance());
            X(i,j) = edit_distance;
            X(j,i) = edit_distance;
            if (num  % 1000000 == 0) {
                cliques::output(graph.id(n1), graph.id(n2), edit_distance, num, ":",total);
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

    int i=0;
    int j=0;
    int num = 0;
    for (auto itr1 = partitions.begin(); itr1 != partitions.end(); ++itr1) {
    	j = i+1;
        for (auto itr2 = itr1; ++itr2 != partitions.end();) {
        	auto p1 = *itr1;
            auto p2 = *itr2;
            cliques::Hungarian hungarian(p1,p2);
            float edit_distance = float(hungarian.edit_distance());
            X(i,j) = edit_distance;
            X(j,i) = edit_distance;
            if (num  % 1000000 == 0) {
                cliques::output(i,j, edit_distance, num, ":", total);
            }
            num++;
            j++;
        }i++;
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

    int i=0;
    for (auto l_itr = landmarks.begin(); l_itr != landmarks.end(); ++l_itr) {
        int j=0;
        for (auto n_itr = landmarks.begin(); n_itr != landmarks.end(); ++n_itr) {
            auto p1 = *n_itr;
            auto p2 = *l_itr;
            cliques::Hungarian hungarian(p1,p2);
            float edit_distance = float(hungarian.edit_distance());
            D_l(i,j) = edit_distance;
            ++j;
        }
        ++i;
    }
    return D_l;
}

/**
 @brief  Embed pairwise distance matrix into another (lower) dimension
 //TODO deal with negative eigen values
 */
arma::mat embed_mds(arma::mat &X,int num_dimen) {
    int N = X.n_cols;
    arma::vec means = arma::zeros<arma::vec>(N);
    for (int i=0; i<N;++i) {
        means(i) = arma::mean(X.col(i));
    }
    double mean = arma::mean(means);
    arma::mat B_n(N,N);
    B_n.zeros();

    for (int i=0;i<N-1;++i) {
        for (int j=i+1;j<N;++j) {
            B_n(i,j) = -(X(i,j) - means(i) - means(j) +mean)/2.0;
            B_n(j,i) = B_n(i,j);
        }
    }

    arma::vec eigvals;
    arma::mat eigvecs;
    eig_sym(eigvals, eigvecs, B_n);

    arma::uvec indices = arma::sort_index(eigvals, 1);
    arma::mat L(num_dimen, N);
    for (int i=0;i<num_dimen;++i) {
        for (int j=0;j<N;++j) {
            L(i,j) = eigvecs(j, indices(i)) * std::sqrt(eigvals(indices(i)));
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
    for (int i=0; i<N;++i) {
        means(i) = arma::mean(D_n.col(i));
    }
    double mean = arma::mean(means);
    arma::mat B_n(N,N);
    B_n.zeros();

    for (int i=0;i<N-1;++i) {
        for (int j=i+1;j<N;++j) {
            B_n(i,j) = -(D_n(i,j) - means(i) - means(j) +mean)/2.0;
            B_n(j,i) = B_n(i,j);
        }
    }

    arma::vec eigvals;
    arma::mat eigvecs;
    eig_sym(eigvals, eigvecs, B_n);

    arma::uvec indices = arma::sort_index(eigvals, 1);
    arma::mat L(num_dimen, N);
    for (int i=0;i<num_dimen;++i) {
        for (int j=0;j<N;++j) {
            L(i,j) = eigvecs(j, indices(i)) * std::sqrt(eigvals(indices(i)));
        }
    }

    arma::mat L_sharp(num_dimen, N);
    for (int i=0;i<num_dimen;++i) {
        for (int j=0;j<N;++j) {
            L_sharp(i,j) = eigvecs(j, indices(i)) / std::sqrt(eigvals(indices(i)));
        }
    }

    int N_all = D_l.n_cols;
    arma::mat X(num_dimen, N_all);
    for (int i=0;i<N_all;++i) {
        arma::vec second_term = D_l.col(i) - means;
        X.col(i) = -L_sharp * second_term/2.0;
    }

    return X;
}

/**
 @brief  Embeds a graph in a lower dim space (convenience function)
 */
template <typename G, typename M>
arma::mat embed_graph(G& graph, M& weights, int num_dim) {
    typedef typename G::Node Node;
    typedef typename G::NodeIt NodeIt;

    // Use all nodes
    std::vector<Node> landmark_nodes;
    for (NodeIt n(graph); n!= lemon::INVALID; ++n) {
        Node node = n;
        landmark_nodes.push_back(node);
    }

    cliques::output("choosing nodes");
    auto X = cliques::find_geodesic_dists(graph, landmark_nodes, weights);
    X.print("X");

    cliques::output("finding embedding");
    auto L = cliques::embed_mds(X, num_dim);
    L.print("L");

    cliques::output("saving");
    arma::mat L_t = arma::trans(L);

    auto D_y = cliques::euclid_pairwise_dists(L_t);
    cliques::output("residual variance", cliques::residual_variance(X, D_y));

    return L_t;
}

template <typename G, typename M, typename T>
void dists_to_graph(G &graph, M &weights, arma::mat &dists, T threshold) {
    typedef typename G::Node Node;
    typedef typename G::Edge Edge;
    int N = dists.n_cols;
    for (int i=0; i<N; ++i) {
        graph.addNode();
    }

    for (int i=0;i<N-1;++i) {
        for (int j=i+1;j<N;++j) {
            if (dists(i,j) == threshold) {
                Node u = graph.nodeFromId(i);
                Node v = graph.nodeFromId(j);
                Edge new_edge = graph.addEdge(u,v);
                weights.set(new_edge, dists(i,j));
            }
        }
    }
}

}
