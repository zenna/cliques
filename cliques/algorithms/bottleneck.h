/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010- 2011 */
#pragma once

#include <lemon/connectivity.h>
#include <map>
#include <set>
#include <lemon/maps.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/smart_graph.h>

#include <cliques/structures/common.h>

namespace clq {

/**
 @brief  Finds the Bottleneck Shortest Path between two points

 @param[in] g The input graph
 @param[in] weights the weights of edges
 @param[in] s the source node
 @param[in] filter this is an edge filter map, required due to recursion.  It should be set such that all nods are valid
 */
float find_bottleneck(
        lemon::SmartGraph &g,
        lemon::IterableValueMap<lemon::SmartGraph, lemon::SmartGraph::Edge,
                float> &weights, lemon::SmartGraph::Node &s,
        lemon::SmartGraph::Node &t, lemon::SmartGraph::EdgeMap<bool> &filter) {

    // Check if theres a single edge left
    lemon::FilterEdges<lemon::SmartGraph> subgraph(g, filter);
    int num_edges = lemon::countEdges(subgraph);
    //int num_edges = countEdges(g);
    if (num_edges == 1) {
        for (lemon::IterableValueMap<lemon::SmartGraph,
                lemon::SmartGraph::Edge, float>::ValueIt itr =
                weights.beginValue(); itr != weights.endValue(); ++itr) {
            for (lemon::IterableValueMap<lemon::SmartGraph,
                    lemon::SmartGraph::Edge, float>::ItemIt item_itr(weights,
                    *itr); item_itr != lemon::INVALID; ++item_itr) {
                if (filter[item_itr] != false) {
                    return *itr;
                }
            }
        }
    }

    std::vector<int> deleted_edges;
    std::vector<lemon::SmartGraph::Node> deleted_u;
    std::vector<lemon::SmartGraph::Node> deleted_v;
    std::vector<float> deleted_weights;
    int num_deleted = 0;

    //cout << "Edges before pruning: " << num_edges << endl;
    //IterableValueMap<filterEdges, Edge, float> newwiehgts(subgraph);
    //FilterEdges<SmartGraph>::EdgeIt a(subgraph);

    for (lemon::IterableValueMap<lemon::SmartGraph, lemon::SmartGraph::Edge,
            float>::ValueIt itr = weights.beginValue(); itr
            != weights.endValue(); ++itr) {
        for (lemon::IterableValueMap<lemon::SmartGraph,
                lemon::SmartGraph::Edge, float>::ItemIt item_itr(weights, *itr); item_itr
                != lemon::INVALID; ++item_itr) {
            if (num_deleted >= num_edges / 2) {
                break;
            }
            if (filter[item_itr] != false) {
                filter[item_itr] = false;
                deleted_edges.push_back(g.id(item_itr));
                num_deleted++;
                deleted_u.push_back(g.u(item_itr));
                deleted_v.push_back(g.v(item_itr));
                deleted_weights.push_back(*itr);
            }
        }
    }
    //for (int i = 0;i<deleted_edges.size();++i) {
    //    g.erase(g.edgeFromId(deleted_edges[i]));
    //}

    //num_edges = countEdges(subgraph);
    num_edges = lemon::countEdges(g);
    //cout << "Edges after pruning: " << num_edges << endl;
    //SubGraph::NodeMap<int> comp
    lemon::SmartGraph::NodeMap<int> comp(g);
    //int num_components = connectedComponents(g,comp);
    int num_components = lemon::connectedComponents(subgraph, comp);

    if (comp[s] != comp[t]) {
        //cout << "ST disjoint - num components = " << num_components << endl;
        lemon::SmartGraph temp_graph;
        for (int i = 0; i < num_components; ++i) {
            temp_graph.addNode();
        }
        s = temp_graph.nodeFromId(comp[s]);
        t = temp_graph.nodeFromId(comp[t]);

        //cout << "Creating new map";
        //std::map<Edge,float> *temp_weight_map = new std::map<Edge,float>;
        clq::setumap edges_map;
        for (unsigned int i = 0; i < deleted_v.size(); ++i) {
            int component_of_v = comp[deleted_v[i]];
            int component_of_u = comp[deleted_u[i]];

            if (component_of_u != component_of_v) {
                std::set<int> components;
                components.insert(component_of_u);
                components.insert(component_of_v);

                clq::setumap::iterator itr = edges_map.find(components);
                if (itr == edges_map.end()) {
                    edges_map[components] = deleted_weights[i];
                } else {
                    if (itr->second < deleted_weights[i]) {
                        itr->second = deleted_weights[i];
                    }
                }
            }
        }
        std::vector<lemon::SmartGraph::Edge> edges;
        for (clq::setumap::iterator itr = edges_map.begin(); itr
                != edges_map.end(); ++itr) {
            std::vector<int> uv;
            for (std::set<int>::iterator sitr = itr->first.begin(); sitr
                    != itr->first.end(); ++sitr) {
                uv.push_back(*sitr);
            }
            edges.push_back(
                    temp_graph.addEdge(temp_graph.nodeFromId(uv[0]),
                            temp_graph.nodeFromId(uv[1])));
        }
        lemon::IterableValueMap<lemon::SmartGraph, lemon::SmartGraph::Edge,
                float> weights(temp_graph); // Map must be made after edges are added
        int i = 0;
        for (clq::setumap::iterator itr = edges_map.begin(); itr
                != edges_map.end(); ++itr) {
            weights.set(edges[i], itr->second);
            i++;
        }
        lemon::SmartGraph::EdgeMap<bool> temp_filter_map(temp_graph, true);
        //istGraph::EdgeMap<bool> filter_map;
        //cout << "recurssinn" << endl;
        return find_bottleneck(temp_graph, weights, s, t, temp_filter_map);
    } else {
        //cout << "Recursing with same map" << endl;
        //cout << "recursing oldskool" << endl;
        return find_bottleneck(g, weights, s, t, filter);
    }
}

}
