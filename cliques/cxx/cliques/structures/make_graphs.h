/* Copyright (c) Modified by Zenna Tavares-zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_MAKEGRAPHS_H
#define CLIQUES_MAKEGRAPHS_H
#include <assert.h>
#include <cstdlib>
#include <iostream>

namespace cliques {

template<typename G>
void make_complete_graph(G &graph, int num_nodes) {
    typedef typename G::NodeIt NodeIt;

    for (int i = 0; i < num_nodes; ++i) {
        graph.addNode();
    }

    for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
        for (NodeIt n2 = n1; ++n2 != lemon::INVALID;) {
            graph.addEdge(n1,n2);
        }
    }
}

template<typename G>
void make_path_graph(G &graph, int num_nodes) {
    typedef typename G::NodeIt NodeIt;

    for (int i = 0; i < num_nodes; ++i) {
        graph.addNode();
    }
    for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
        NodeIt n2 = n1;
        ++n2;
        if (n2 == lemon::INVALID) {
            break;
        }
        graph.addEdge(n1,n2);
    }
}

}
#endif