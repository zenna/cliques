#pragma once

#include <iostream>
#include <vector>
#include <cmath>


extern "C" void dgpadm_(int* ideg, int* m, double* t, double* A, int* ldh,
        double* wsp, int* lwsp, int* iwsp, int* iexp, int* ns, int* iflag);

namespace clq {

std::vector<double> exp(std::vector<double> matrix, double t, int order) {
    int ideg = 6;
    int m = order;

    int lda = order;
    int ldh = lda;
    int lwsp = 4 * ldh * ldh + ideg + 1;
    double *wsp = new double[lwsp];
    double *A = &matrix.front();
    int *iwsp = new int[ldh];


    // output
    int iexp, ns, iflag;
    //TODO: This does not work with MATLAB! But no idea why!
    dgpadm_(&ideg, &m, &t, A, &lda, wsp, &lwsp, iwsp, &iexp, &ns, &iflag);

    double *start = wsp + iexp -1;

    std::vector<double> output;
    for (unsigned int i = 0; i < matrix.size(); ++i) {
        double a = start[i];
        output.push_back(a);
    }

    delete[] iwsp;
    delete[] wsp;

    return output;
}

std::vector<double> create_exponential_markov_times(double start_time,
        double end_time, int num_steps) {
    std::vector<double> markov_times;
    double start_log = std::log(start_time);
    double end_log = std::log(end_time);

    double increment = (end_log - start_log) / float(num_steps - 1);
    for (int i = 0; i < num_steps; ++i) {
        double current_log = start_log + i * increment;
        markov_times.push_back(std::exp(current_log));
    }
    return markov_times;
}

double discrete_gauss_kernel(int N, double T) {
    return 1.0 / (N + 1);
}

template<typename G, typename M>
void graph_to_exponential_graph(G &graph, M &weights, G &exp_graph, M &exp_graph_weights, double markov_time) {
    std::vector<double> minus_t_D_inv_L;
    std::vector<double> node_weighted_degree;

    typedef typename G::Node Node;
    typedef typename G::NodeIt NodeIt;
    typedef typename G::EdgeIt EdgeIt;
    /////////////////
    // read out graph into matrix
    //////////////////

    // number of nodes N
    int N = lemon::countNodes(graph);
    // -D^-1*L*t == t(B-I)
    minus_t_D_inv_L = std::vector<double>(N * N, 0);
    node_weighted_degree = std::vector<double>(N, 0);

    // get weighted degree of nodes
    for (int i = 0; i < N; ++i) {
        typename G::Node temp_node = graph.nodeFromId(i);
        node_weighted_degree[i] = find_weighted_degree(graph, weights,
                temp_node);
    }

    //initialise matrix and set diagonals to minus identity
    for (int i = 0; i < N * N; ++i) {
        if (i % (N + 1) == 0) {
            minus_t_D_inv_L[i] = -1;
        } else {
            minus_t_D_inv_L[i] = 0;
        }

    }

    // fill in the rest
    for (EdgeIt e(graph); e != lemon::INVALID; ++e) {
        Node u = graph.u(e);
        Node v = graph.v(e);
        int node_id_u = graph.id(u);
        int node_id_v = graph.id(v);
        double weight_uv = weights[e];
        // B_uv
        minus_t_D_inv_L[node_id_v + N * node_id_u] += weight_uv
                / node_weighted_degree[node_id_v];
        // B_vu
        minus_t_D_inv_L[node_id_u + N * node_id_v] += weight_uv
                / node_weighted_degree[node_id_u];
    }

    // call expokit
    int N2 = node_weighted_degree.size();
    std::vector<double> exp_graph_vec = clq::exp(minus_t_D_inv_L,
            markov_time, N2);

    // reserve memory space for number of nodes
    exp_graph.reserveNode(N2);
    exp_graph.reserveEdge(N2 + (N2 * (N2 - 1)) / 2);

    // add nodes
    for (int i = 0; i < N2; ++i) {
        exp_graph.addNode();
    }

    for (int i = 0; i < N2; ++i) {
        for (int j = i; j < N2; ++j) {
            double weight = exp_graph_vec[N2 * i + j] * node_weighted_degree[j];
            if (weight > 0) {
                lemon::SmartGraph::Edge edge = exp_graph.addEdge(
                        exp_graph.nodeFromId(i), exp_graph.nodeFromId(j));
                exp_graph_weights.set(edge, weight);
            }
        }
    }
}

}
