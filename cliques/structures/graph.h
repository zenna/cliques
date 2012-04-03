#ifndef CLIQUES_GRAPH_H
#define CLIQUES_GRAPH_H

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

namespace cliques {
    template <class T>
    struct Graph
    {
        T graph;

        float A(int node1_id, int node2_id) {
            lemon::ListGraph::Node n1 = this->graph.nodeFromId(node1_id);
            lemon::ListGraph::Node n2 = this->graph.nodeFromId(node2_id);

            for(lemon::ListGraph::IncEdgeIt e(this->graph, n1); e!= lemon::INVALID; ++e) {
                //std::cout << "node1 " << this->graph.id(n1) << " node2 " << this->graph.id(n2) << std::endl;
                //std::cout << "base: " << this->graph.id(graph.baseNode(e)) <<  " running: " << this->graph.id(graph.runningNode(e)) << std::endl;
                if (this->graph.runningNode(e) == n2) {
                    //std::cout << "one" << std::endl;
                    return 1.0;
                }
            }
            //std::cout << "hero" << std::endl;
            return 0.0;
        };

        int k(int node_id) {
            int count = 0;
            lemon::ListGraph::Node n = graph.nodeFromId(node_id);

            for(lemon::ListGraph::IncEdgeIt e(graph, n); e!= lemon::INVALID; ++e) {
                ++count;
            }
            return count;
        };

        int num_edges() {return lemon::countEdges(this->graph);}
        int num_nodes() {return lemon::countNodes(this->graph);}
    };
} //namespace cliques

#endif //CLIQUES_STRUCTURES_H
