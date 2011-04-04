#include <map>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/random.h>
#include <lemon/lgf_writer.h>

#include <lemon/concepts/digraph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/smart_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/static_graph.h>

#include <vector>


#if defined USE_BOOST
#include <boost/unordered_map.hpp>
#else
#include <unordered_map>
#endif

using namespace lemon;
using namespace std;

typedef SmartGraph::Node Node;
typedef SmartGraph::Edge Edge;
typedef SmartGraph::EdgeMap<int> LengthMap;
typedef lemon::SmartGraph Graph;


struct set_hash
{
    size_t operator()(std::set<int> set_key) const
    {
        int alpha = 0;
        for (set<int>::iterator itr = set_key.begin(); itr != set_key.end(); ++itr)
          alpha += *itr;
        return alpha; //shash()(bitset.to_string());
    }
};

#if defined USE_BOOST
typedef boost::unordered_map<std::set<int>,float,set_hash> umap;
#else
typedef std::unordered_map<std::set<int>,float,set_hash> umap;
#endif


//Speedups:
//Remove unnecessary compuatiaton
//Use subgraph instead of graph copy

void print_map(std::map<int,int> my_map, bool add_one = false) {

    for (std::map<int,int>::iterator itr = my_map.begin(); itr != my_map.end(); ++itr) {
        int node = itr->first;
        int part = itr->second;
        if (add_one) {
            node++;
            part++;
        }
        std::cout << node << ":" << part << " ";
    }
    std::cout << "\n";
}

////template <class float>
//float find_bottleneck (ListGraph &g, IterableValueMap<ListGraph, Edge, float> &weights, Node &s, Node &t) {
//    // Check if theres a single edge left
//    int num_edges = countEdges(g);
//    if (num_edges == 1) {
//        IterableValueMap<ListGraph, Edge, float>::ValueIt itr = weights.beginValue();
//        return *itr;
//    }
//
//    vector<int> deleted_edges;
//    vector<Node> deleted_u;
//    vector<Node> deleted_v;
//    vector<float> deleted_weights;
//    int num_deleted = 0;
//
//    cout << "Edges before pruning: " << num_edges << endl;
//    for (IterableValueMap<ListGraph, Edge, float>::ValueIt itr = weights.beginValue();itr != weights.endValue(); ++itr) {
//        for (IterableValueMap<ListGraph, Edge, float>::ItemIt item_itr(weights, *itr); item_itr != INVALID; ++item_itr) {
//            if (num_deleted >= num_edges / 2) {
//                break;
//            }
//            deleted_edges.push_back(g.id(item_itr));
//            num_deleted++;
//            deleted_u.push_back(g.u(item_itr));
//            deleted_v.push_back(g.v(item_itr));
//            deleted_weights.push_back(*itr);
//
//        }
//    }
//
//    for (int i = 0;i<deleted_edges.size();++i) {
//        g.erase(g.edgeFromId(deleted_edges[i]));
//    }
//
//    num_edges = countEdges(g);
//    cout << "Edges after pruning: " << num_edges << endl;
//    Graph::NodeMap<int> comp(g);
//    int num_components = connectedComponents(g,comp);
//
//    if (comp[s] != comp[t]) {
//        cout << "ST disjoint - num components = " << num_components << endl;
//        ListGraph temp_graph;
//        for (int i=0;i<num_components;++i) {
//            temp_graph.addNode();
//        }
//        s = temp_graph.nodeFromId(comp[s]);
//        t = temp_graph.nodeFromId(comp[t]);
//
//        cout << "Creating new map";
//        std::map<Edge,float> *temp_weight_map = new std::map<Edge,float>;
//        for (unsigned int i = 0;i<deleted_v.size();++i) {
//            if (i%100 == 0) {
//                cout << "Checking edge " << i << " out of "<< deleted_v.size() << endl;
//            }
//            int component_of_v = comp[deleted_v[i]];
//            int component_of_u = comp[deleted_u[i]];
//            if (component_of_u != component_of_v) {
//                Node node_u = temp_graph.nodeFromId(component_of_u);
//                Node node_v = temp_graph.nodeFromId(component_of_v);
//                // Check if there is an edge
//                Edge e = findEdge(temp_graph,node_u,node_v);
//                if (e != INVALID) {
//                    if ((*temp_weight_map)[e] < deleted_weights[i]) {
//                        (*temp_weight_map)[e] == deleted_weights[i];
//                    }
//                }
//                else {
//                    Edge e = temp_graph.addEdge(node_u,node_v);
//                    (*temp_weight_map)[e] = deleted_weights[i];
//                }
//            }
//        }
//        cout << "Checking copying stdmap to weight map" << endl;
//        IterableValueMap<ListGraph, Edge, float> weights(temp_graph);
//        for (std::map<Edge,float>::iterator itr = temp_weight_map->begin(); itr != temp_weight_map->end(); ++itr) {
//            weights.set(itr->first,itr->second);
//        }
//        delete temp_weight_map;
//        cout << "recursing with new CC map" << endl;
//        return find_bottleneck(temp_graph,weights,s,t);
//    }
//    else {
//        cout << "Recursing with same map" << endl;
//        return find_bottleneck(g,weights,s,t);
//    }
//}


float find_bottleneck (SmartGraph &g, IterableValueMap<SmartGraph, Edge, float> &weights, Node &s, Node &t, SmartGraph::EdgeMap<bool> &filter) {
    // Check if theres a single edge left
    FilterEdges<SmartGraph> subgraph(g, filter);
    int num_edges = countEdges(subgraph);
    //int num_edges = countEdges(g);
    if (num_edges == 1) {
        for (IterableValueMap<SmartGraph, Edge, float>::ValueIt itr = weights.beginValue();itr != weights.endValue(); ++itr) {
            for (IterableValueMap<SmartGraph, Edge, float>::ItemIt item_itr(weights, *itr); item_itr != INVALID; ++item_itr) {
                if (filter[item_itr] != false) {
                    return *itr;
                }
            }
        }
    }

    cout << "num_nodes = " << countNodes(subgraph) << endl;

    vector<int> deleted_edges;
    vector<Node> deleted_u;
    vector<Node> deleted_v;
    vector<float> deleted_weights;
    int num_deleted = 0;

    cout << "Edges before pruning: " << num_edges << endl;
    //IterableValueMap<filterEdges, Edge, float> newwiehgts(subgraph);
    //FilterEdges<SmartGraph>::EdgeIt a(subgraph);

    for (IterableValueMap<SmartGraph, Edge, float>::ValueIt itr = weights.beginValue();itr != weights.endValue(); ++itr) {
        for (IterableValueMap<SmartGraph, Edge, float>::ItemIt item_itr(weights, *itr); item_itr != INVALID; ++item_itr) {
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

    num_edges = countEdges(subgraph);
    //num_edges = countEdges(g);
    cout << "Edges after pruning: " << num_edges << endl;
    //SubGraph::NodeMap<int> comp
    Graph::NodeMap<int> comp(g);

    for (IterableValueMap<SmartGraph, Edge, float>::ValueIt itr = weights.beginValue();itr != weights.endValue(); ++itr) {
        for (IterableValueMap<SmartGraph, Edge, float>::ItemIt item_itr(weights, *itr); item_itr != INVALID; ++item_itr) {
            if (filter[item_itr] != false) {
                cout << g.id(g.u(item_itr)) << "-" << g.id(g.v(item_itr)) << ":" << *itr << endl;
            }
        }
    }
    //int num_components = connectedComponents(g,comp);
    int num_components = connectedComponents(subgraph,comp);

    if (comp[s] != comp[t]) {
        cout << "ST disjoint - num components = " << num_components << endl;
        SmartGraph temp_graph;
        for (int i=0;i<num_components;++i) {
            temp_graph.addNode();
        }
        s = temp_graph.nodeFromId(comp[s]);
        t = temp_graph.nodeFromId(comp[t]);

        //cout << "Creating new map";
        //std::map<Edge,float> *temp_weight_map = new std::map<Edge,float>;
        umap edges_map;
        for (unsigned int i = 0;i<deleted_v.size();++i) {
            int component_of_v = comp[deleted_v[i]];
            int component_of_u = comp[deleted_u[i]];

            if (component_of_u != component_of_v) {
                std::set<int> components;
                components.insert(component_of_u);
                components.insert(component_of_v);

                umap::iterator itr = edges_map.find(components);
                if (itr == edges_map.end()) {
                     edges_map[components] = deleted_weights[i];
                }
                else {
                    if (itr->second < deleted_weights[i]) {
                        itr->second = deleted_weights[i];
                    }
                }
            }
        }
        vector<Edge> edges;
        for (umap::iterator itr = edges_map.begin(); itr != edges_map.end(); ++itr) {
            std::vector<int> uv;
            for (std::set<int>::iterator sitr = itr->first.begin(); sitr != itr->first.end(); ++sitr) {
                uv.push_back(*sitr);
            }
            edges.push_back(temp_graph.addEdge(temp_graph.nodeFromId(uv[0]),temp_graph.nodeFromId(uv[1])));
        }
        IterableValueMap<SmartGraph, Edge, float> weights(temp_graph); // Map must be made after edges are added
        int i = 0;
        for (umap::iterator itr = edges_map.begin(); itr != edges_map.end(); ++itr) {
            weights.set(edges[i],itr->second);
            i++;
        }
        SmartGraph::EdgeMap<bool> temp_filter_map(temp_graph, true);
        //istGraph::EdgeMap<bool> filter_map;
        //cout << "recurssinn" << endl;
        return find_bottleneck(temp_graph,weights,s,t,temp_filter_map);
    }
    else {
        cout << "Recursing with same map" << endl;
        //cout << "recurssinn oldskool" << endl;
        return find_bottleneck(g,weights,s,t,filter);
    }
}

// Save as set of bitsets/Set of maps/
// Key 0 1 1 0 2 0
// Key Map
#include <string>
#include <map>

#include <iostream>
#include <fstream>
void save_partition(std::map<int, std::map<int, int> > partitions_map, std::string filename) {
    ofstream partition_file;
    partition_file.open(filename.c_str());
    for (std::map<int, std::map<int, int> >::iterator itr = partitions_map.begin(); itr != partitions_map.end(); ++itr) {
        partition_file << itr->first  << ",";
        for (std::map<int, int>::iterator part_itr = itr->second.begin(); part_itr != itr->second.end(); ++part_itr) {
            partition_file << part_itr->first << "," << part_itr->second << ",";
        }
        partition_file << std::endl;
    }
    partition_file.close();
}


int main()
{
    char * UNIPROJ_HOME;
    UNIPROJ_HOME= getenv ("UNPROJ_HOME");
    //if (pPath!=NULL)
    // Load neighbour graph
    // Load stabilites into array 2MB * # Samples - parameterise to do different range
    // For each stability
        // Find maxima using arrayt as node value
        // Generate edge map from in nodes
        // BSP to find  pairwise stabilities -> store to file

    SmartGraph g;
    for (int i =0; i<12; ++i) {
        g.addNode();
    }
    std::vector<float> stabilities;
    stabilities.push_back(100);
    stabilities.push_back(30);
    stabilities.push_back(90);
    stabilities.push_back(59);
    stabilities.push_back(12);
    stabilities.push_back(80);
    stabilities.push_back(44);
    stabilities.push_back(58);
    stabilities.push_back(64);
    stabilities.push_back(60);
    stabilities.push_back(100);
    stabilities.push_back(80);

    g.addEdge(g.nodeFromId(0),g.nodeFromId(1));
    g.addEdge(g.nodeFromId(1),g.nodeFromId(2));
    g.addEdge(g.nodeFromId(0),g.nodeFromId(2));
    g.addEdge(g.nodeFromId(1),g.nodeFromId(10));
    g.addEdge(g.nodeFromId(2),g.nodeFromId(3));
    g.addEdge(g.nodeFromId(3),g.nodeFromId(4));
    g.addEdge(g.nodeFromId(3),g.nodeFromId(5));
    g.addEdge(g.nodeFromId(4),g.nodeFromId(5));
    g.addEdge(g.nodeFromId(6),g.nodeFromId(7));
    g.addEdge(g.nodeFromId(6),g.nodeFromId(8));
    g.addEdge(g.nodeFromId(7),g.nodeFromId(8));
    g.addEdge(g.nodeFromId(4),g.nodeFromId(7));
    g.addEdge(g.nodeFromId(5),g.nodeFromId(9));
    g.addEdge(g.nodeFromId(9),g.nodeFromId(11));
    g.addEdge(g.nodeFromId(9),g.nodeFromId(10));
    g.addEdge(g.nodeFromId(10),g.nodeFromId(11));


    IterableValueMap<SmartGraph, Edge, float> weights(g);

    std::map<std::pair<int,int>,float> pairwise_bottlenecks;
    for (SmartGraph::EdgeIt itr(g); itr != INVALID; ++itr) {
        float u_stab = stabilities[g.id(g.u(itr))];
        float v_stab = stabilities[g.id(g.v(itr))];
        float lesser_stab = u_stab < v_stab ? u_stab : v_stab;
        weights.set(itr,lesser_stab);
    }

    for (IterableValueMap<SmartGraph, Edge, float>::ValueIt itr = weights.beginValue();itr != weights.endValue(); ++itr) {
        for (IterableValueMap<SmartGraph, Edge, float>::ItemIt item_itr(weights, *itr); item_itr != INVALID; ++item_itr) {
            cout << g.id(g.u(item_itr)) << "-" << g.id(g.v(item_itr)) << ":" << *itr << endl;
        }
    }

    SmartGraph::EdgeMap<bool> filter(g, true);

    Node source = g.nodeFromId(0);
    Node target = g.nodeFromId(6);

    cout << find_bottleneck (g, weights, source , target, filter);

    /*Node v0=g.addNode();
    Node v1=g.addNode();
    Node v2=g.addNode();
    Node v3=g.addNode();
    Node v4=g.addNode();
    Node v5=g.addNode();
    Node v6=g.addNode();

    Edge v0_v1=g.addEdge(v0,v1);
    Edge v1_v2=g.addEdge(v1,v2);
    Edge v0_v3=g.addEdge(v0, v3);
    Edge v3_v2=g.addEdge(v3, v2);
    Edge v3_v1=g.addEdge(v3, v1);
    Edge v0_v4=g.addEdge(v0,v4);
    Edge v4_v2=g.addEdge(v4,v2);
    Edge v0_v5=g.addEdge(v0,v5);
    Edge v0_v6=g.addEdge(v0,v6);
    Edge v6_v2=g.addEdge(v6,v2);
    Edge v5_v6=g.addEdge(v5,v6);

    IterableValueMap<ListGraph, Edge, float> cont_weights(g);

    cont_weights.set(v0_v1,-0.3);
    cont_weights.set(v1_v2,0.9);
    cont_weights.set(v0_v3,0.5);
    cont_weights.set(v3_v2,1.0);
    cont_weights.set(v3_v1,0.55);
    cont_weights.set(v0_v4,-0.2);
    cont_weights.set(v4_v2,-0.3);
    cont_weights.set(v0_v5,0.6);
    cont_weights.set(v0_v6,.4);
    cont_weights.set(v6_v2,0.8);
    cont_weights.set(v5_v6,0.6);

    ListGraph rg;

    cout << "Generating Nodes" << endl;

    int size = 10000;
    for (int i =0;i<size;++i) {
        rg.addNode();
    }

    getchar();

    cout << "Generating Edges" << endl;

    for (ListGraph::NodeIt itr(rg); itr != INVALID; ++itr) {
        for (int i=0;i<32;++i) {
            int v = rnd[size];
            Node nv = g.nodeFromId(v);
            if (itr != nv)
                rg.addEdge(itr, nv);
        }
    }
    getchar();

    graphWriter(rg, "graph2.lgf").run();

    cout << "Generating Weights " << endl;

    IterableValueMap<ListGraph, Edge, float> weights(rg);
    for (ListGraph::EdgeIt itr(rg); itr != INVALID; ++itr) {
        float a = rnd();
        weights.set(itr,a);
    }

    Node source = rg.nodeFromId(5);
    Node target = rg.nodeFromId(size/2);
    getchar();

    cout << "Connected? : ";
    cout << connected(rg) << endl;

    cout << "finding bneck";

    cout << "bottleneck is: " << find_bottleneck(rg,weights,source,target) << endl;

    // Save Neighbours
    graphWriter(rg, "graph.lgf").run();


    // Save Partitions
    // Save Modularity/Stability

    //Program
    /*load graph from file
    load quality scores into map key: vector or array (2 mb)

    for number of quality measures
        follow neighbours to find maxima using values stored in array
        generate graph weights by choosing min of thingy
        find bottleneck path
        save output*/
    return 0;
}
