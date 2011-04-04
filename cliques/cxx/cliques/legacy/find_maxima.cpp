#include <map>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/concepts/digraph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/smart_graph.h>
#include <lemon/random.h>
#include <lemon/lgf_writer.h>
#include <lemon/lgf_reader.h>
#include <lemon/static_graph.h>

#include <vector>
#include <set>
#include <cstdio>

#include <lemon/random.h>

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

void reorder_ints(int &i, int &j) {
    int buffer;
    if (i <= j) {
        return;
    }
    else {
        buffer = i;
        i = j;
        j = buffer;
    }
}

void reorder_map(std::map<int,int> &partition_map) {
    int start_part = 0;
    std::map<int,int> temp_partition_map = partition_map;

    for (std::map<int,int>::iterator itr = temp_partition_map.begin(); itr != temp_partition_map.end(); ++itr) {
        int reference_part = itr->second;
	    std::vector<int> to_remove;

        for (std::map<int,int>::iterator new_itr = temp_partition_map.begin(); new_itr != temp_partition_map.end(); ++new_itr) {
            if (new_itr->second == reference_part) {
                partition_map[new_itr->first] = start_part;

                if (itr->first != new_itr->first)
                   to_remove.push_back(new_itr->first);
            }
        }
        for (unsigned int i =0; i< to_remove.size(); ++i) {
            temp_partition_map.erase(to_remove[i]);
        }
        start_part++;
    }
}


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

    vector<int> deleted_edges;
    vector<Node> deleted_u;
    vector<Node> deleted_v;
    vector<float> deleted_weights;
    int num_deleted = 0;

    //cout << "Edges before pruning: " << num_edges << endl;
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

    //num_edges = countEdges(subgraph);
    num_edges = countEdges(g);
    //cout << "Edges after pruning: " << num_edges << endl;
    //SubGraph::NodeMap<int> comp
    Graph::NodeMap<int> comp(g);
    //int num_components = connectedComponents(g,comp);
    int num_components = connectedComponents(subgraph,comp);

    if (comp[s] != comp[t]) {
        //cout << "ST disjoint - num components = " << num_components << endl;
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
        //cout << "Recursing with same map" << endl;
        //cout << "recurssinn oldskool" << endl;
        return find_bottleneck(g,weights,s,t,filter);
    }
}

void save_vector(std::vector<float> vec_to_save, std::string filename) {
    ofstream vector_file;
    vector_file.open(filename.c_str());
    //partition_file << partition_map.size() << endl;
    for (std::vector<float>::iterator itr = vec_to_save.begin(); itr != vec_to_save.end(); ++ itr) {
        vector_file << *itr << endl;
    }
}


// Save as set of bitsets/Set of maps/
// Key 0 1 1 0 2 0
// Key Map
#include <string>
#include <map>

#include <iostream>
#include <fstream>
void save_maxima(std::vector<std::set<int> > &all_maxima, std::string filename,std::string num_filename) {
    ofstream maxima_num_file;
    maxima_num_file.open(num_filename.c_str());
    ofstream maxima_file;
    maxima_file.open(filename.c_str());

    for (std::vector<set<int> >::iterator a_itr = all_maxima.begin(); a_itr != all_maxima.end(); ++a_itr) {
        //maxima_file<< a_itr->size() << ",";
        maxima_num_file << a_itr->size() << endl;
        for (std::set<int>::iterator itr = a_itr->begin(); itr != a_itr->end(); ++itr) {
            maxima_file<< *itr << ",";
        }
        maxima_file << endl;
    }
    maxima_num_file.close();
    maxima_file.close();
}

void load_maxima(std::vector<std::set<int> > &all_maxima, std::string filename) {
    ifstream maxima_file(filename.c_str());
    std::string line;
    std::string maxima;

    if ( !maxima_file.is_open() ) {
        std::cout << "couldn't open file" << endl;
        exit(1);
    }
    while (std::getline(maxima_file,line)) {
        std::stringstream lineStream(line);
        std::set<int> current_maxima;
        while (std::getline(lineStream,maxima,',')) {
            current_maxima.insert(atoi(maxima.c_str()));
        }
        all_maxima.push_back(current_maxima);
    }
    maxima_file.close();
}

void load_stabilities(std::string filename, int start_sample, int end_sample, float **stabilities, int  &num_partitions) {
    ifstream stabilities_file(filename.c_str());
    std::string line;
    std::string stability;

    if ( !stabilities_file.is_open() ) {
        std::cout << "couldn't open file" << endl;
        exit(1);
    }
    std::getline(stabilities_file,line);
    num_partitions = atoi(line.c_str());

    *stabilities = new float[num_partitions * (end_sample - start_sample + 1)];

    int sample_in_array = 0;
    while (std::getline(stabilities_file,line)) {
        std::stringstream lineStream(line);
        int sample_num = 0;
        std::getline(lineStream,stability,',');

        int offset = 0;

        while (std::getline(lineStream,stability,',')) {
            if (sample_num < start_sample) {
                sample_num++;
                continue;
            }
            else if (sample_num > end_sample) {
                sample_num++;
                break;
            }
            else {
                (*stabilities)[sample_in_array + offset*num_partitions] = atof(stability.c_str());
                sample_num++;
                offset++;
            }
        }
        sample_in_array++;
    }
    stabilities_file.close();
}

std::vector<int> random_walk(int num_steps, float eta, SmartGraph &g, float **stabilities, float &eta_star, SmartGraph::NodeMap<int> &degrees, int offset) {
    int num_nodes = countNodes(g);
    SmartGraph::Node start = g.nodeFromId(rnd.integer(num_nodes));
    std::vector<int> symbols;
    for (int i=0;i<num_steps;++i) {
        float old_stab = (*stabilities)[g.id(start) + offset];
        int neighbour_num = rnd.integer(degrees[start]);
        SmartGraph::OutArcIt n_itr(g,start);
        for (int i=0;i<neighbour_num;++i) {
            ++n_itr;
        }
        SmartGraph::Node neighbour = g.runningNode(n_itr);
        float new_stab = (*stabilities)[g.id(neighbour) + offset];

        if (abs(new_stab - old_stab) > eta_star) {
            eta_star = abs(new_stab - old_stab);
        }

        if ((new_stab - old_stab) < -eta) {
            symbols.push_back(-1);
        }
        else if (abs(new_stab - old_stab) <= eta) {
            symbols.push_back(0);
        }
        else if ((new_stab - old_stab) > eta) {
            symbols.push_back(1);
        }
        start = neighbour;
    }
    return symbols;
}

#include <math.h>
float calculate_entropy(std::vector<int> &symbols) {
    float H = 0.0;

    for (int p=-1;p<2;++p) {
        for (int q=-1;q<2;++q) {
            if (p!=q) {
                float P = 0;
                for (int i=0;i<symbols.size()-1 ;++i) {
                    if ((symbols[i] == p) && (symbols[i+1] == q)) {
                        P++;
                    }
                }
                P = P/symbols.size();
                float logP = log(P)/log(6);
                if (P!=0)
                    H = H + P*logP;
            }
        }
    }
    return -H;
}

void save_bottlenecks(std::vector<std::map<std::pair<int,int>,float> > &all_bottlenecks, std::string filename) {
    ofstream bottleneck_file;
    bottleneck_file.open(filename.c_str());

    int i = 0;
    for (std::vector<std::map<std::pair<int,int>,float> >::iterator itr = all_bottlenecks.begin(); itr != all_bottlenecks.end(); ++itr) {
        for (std::map<std::pair<int,int>,float>::iterator mitr = itr->begin(); mitr != itr->end(); ++mitr) {
            bottleneck_file << i << "," << mitr->first.first << "," << mitr->first.second << "," << mitr->second << std::endl;
        }
        i++;
    }
    bottleneck_file.close();
};

void load_bottlenecks(std::vector<std::map<std::pair<int,int>,float> > &all_bottlenecks, std::string filename) {
    ifstream bottleneck_file(filename.c_str());
    std::string line;
    std::string value;

    if ( !bottleneck_file.is_open() ) {
        std::cout << "couldn't open bottleneck file" << endl;
        exit(1);
    }
    std::map<std::pair<int,int>,float> temp_bottleneck;
    all_bottlenecks.push_back(temp_bottleneck);

    while (std::getline(bottleneck_file,line)) {
        std::stringstream lineStream(line);
        std::getline(lineStream,value,',');
        int at_time = atoi(value.c_str());
        int diff = at_time + 1 - all_bottlenecks.size();

        for (int i = 0;i<diff;++i) {
            std::map<std::pair<int,int>,float> temp_bottleneck;
            all_bottlenecks.push_back(temp_bottleneck);
        }

        std::getline(lineStream,value,',');
        int maxima_1 = atoi(value.c_str());
        std::getline(lineStream,value,',');
        int maxima_2 = atoi(value.c_str());
        std::getline(lineStream,value,',');
        float bottleneck = atof(value.c_str());
        all_bottlenecks[at_time][make_pair(maxima_1,maxima_2)] = bottleneck;
    }
    bottleneck_file.close();
}

void load_vector(std::vector<float> &output_vector, std::string filename) {
    ifstream input_file(filename.c_str());
    std::string line;

    if ( !input_file.is_open() ) {
        std::cout << "couldn't open file" << endl;
        exit(1);
    }

    while (std::getline(input_file,line)) {
        output_vector.push_back(atof(line.c_str()));
    }
    input_file.close();
}

//
//output, .max
int main(int argc, char* argv[])
{
    std::string UNIPROJ_HOME = std::getenv ("UNIPROJ_HOME");
    //TODO
    //pass params from commandline
    //Stop revisit if problem
    //UIse adaptors if problem
    //Test range of t e.g .2-3
    //Save BSP to file
    //Compute DG graph from this
    //BUG NUM MAXIMA WRONG STILL
    //BUG ERASE
    int num_partitions;
    std::string basename = argv[1];// || "barbell_n16";
    std::vector<set<int> > all_maxima;
    std::vector<float> taus;
    load_vector(taus,UNIPROJ_HOME+"/data/taus/"+basename+".taus");


    int end, start;
    if (argc > 2) {
        end = atoi(argv[3]);
        start = atoi(argv[2]);
    }
    else {
        start = 0;
        end = taus.size() - 1;
    }
    int range = end - start + 1;



    float *stabilities;
    load_stabilities(UNIPROJ_HOME+"/data/stabilities/"+basename+".stabilities",start,end, &stabilities, num_partitions);

    // Temporary To create Dot graph vis
    {
        ofstream dot_file;
        dot_file.open("graph.gv");
        dot_file << "graph {" << endl
            << "node[width=.25,height=.375,fontsize=0,shape=none,label=\"\"]" << endl
            << "node [shape=none]" << endl
            << "edge[penwidth=0.05,label=\"\",len=4000.0]" << endl
            << "K=10000" << endl
            << "bgcolor=\"black\"" << endl;
        cout << "creating dot graph" << endl;
        SmartGraph gr;
        graphReader(gr, UNIPROJ_HOME+"/data/landscapes/"+basename+"_landscape.lgf").run();

        int offset = 320*num_partitions;

        IterableValueMap<SmartGraph, Edge, float> weights(gr);
        float tau = taus[offset/num_partitions];
        cout << "t =  " << tau << endl;
        std::map<std::pair<int,int>,float> pairwise_bottlenecks;


        float min_weight = 1000;
        float max_weight = -1000;
        int self_connected = 0;
        for (SmartGraph::EdgeIt itr(gr); itr != INVALID; ++itr) {
            if (gr.id(gr.u(itr)) != gr.id(gr.v(itr))) {
                float u_stab = stabilities[gr.id(gr.u(itr)) + offset];
                float v_stab = stabilities[gr.id(gr.v(itr)) + offset];
                float lesser_stab = u_stab < v_stab ? u_stab : v_stab;

                //cout << "self connected:" << self_connected << endl;
                //self_connected++;
                if (lesser_stab < min_weight)
                    min_weight = lesser_stab;

                if (lesser_stab > max_weight)
                    max_weight = lesser_stab;

                weights.set(itr,lesser_stab);
            }
            //dot_file << gr.id(gr.u(itr)) << " -- " << gr.id(gr.v(itr)) << "[color=\"" << abs(lesser_stab) << " 0.700 0.700\"]" << endl;
        }

        float shift = 0.0 - min_weight;
        float range = max_weight - min_weight;

        for (SmartGraph::EdgeIt itr(gr); itr != INVALID; ++itr) {
            if (gr.id(gr.u(itr)) != gr.id(gr.v(itr))) {
                float u_stab = stabilities[gr.id(gr.u(itr)) + offset];
                float v_stab = stabilities[gr.id(gr.v(itr)) + offset];
                float lesser_stab = u_stab < v_stab ? u_stab : v_stab;

                if (lesser_stab < min_weight)
                    min_weight = lesser_stab;

                if (lesser_stab > max_weight)
                    max_weight = lesser_stab;

                weights.set(itr,lesser_stab);
                dot_file << gr.id(gr.u(itr)) << " -- " << gr.id(gr.v(itr)) << "[color=\"" << (lesser_stab+shift)/range<< " 1.0 1.0\"]" << endl;
            }
        }

        cout << "min: " << min_weight << " max: " << max_weight << endl;


        /*for (SmartGraph::NodeIt n(gr); n != INVALID; ++n) {
            for (SmartGraph::OutArcIt a(gr,n); a != INVALID; ++a) {
                dot_file << gr.id(n) << " -- " << gr.id(gr.target(a)) << endl;
            }
        }*/
        dot_file << "}";
        dot_file.close();
    }

    return 0;




    {
        cout << "loading" << endl;
        SmartGraph gr;
        graphReader(gr, UNIPROJ_HOME+"/data/landscapes/"+basename+"_landscape.lgf").run();

        cout << "finding maxima" << endl;
        for (int offset = 0; offset < range*num_partitions; offset = offset + num_partitions) {
            set<int> maxima;
            int num_iterations = 0;
            // Repeat this for every markov time t

            for (SmartGraph::NodeIt n(gr); n != INVALID; ++n) {
                // Check part_itr hasn't already been visited
                if (num_iterations % 100 == 0) {
                    //std::cout << "the number of iterations is" << num_iterations << std::endl;
                    //std::cout << "so far number of maxima is " << maxima.size() << std::endl;
                    //std::cout << "size of visited is " << visited.size() << std::endl;
                    //std::cout << "num skipped " << num_skipped << std::endl;
                }
                SmartGraph::Node best_neighbour = n;
                while (1) {
                    bool has_improved = false;
                    float best_score = stabilities[gr.id(best_neighbour) + offset];
                    for (SmartGraph::OutArcIt a(gr,best_neighbour); a != INVALID; ++a) {
                        if (best_score < stabilities[gr.id(gr.target(a)) + offset]) {
                            best_score = stabilities[gr.id(gr.target(a)) + offset];
                            best_neighbour = gr.target(a);
                            has_improved = true;
                        }
                    }


                    if (has_improved == false) {
                        maxima.insert(gr.id(best_neighbour));
                        break;
                    }
                }
                num_iterations++;
            }
            std::cout << maxima.size() << std::endl;
            all_maxima.push_back(maxima);
        }
    }
    save_maxima (all_maxima, UNIPROJ_HOME+"/data/maxima/"+basename+".max",UNIPROJ_HOME+"/data/maxima/"+basename+".nummax");
    //load_maxima(all_maxima,UNIPROJ_HOME+"/data/maxima/"+basename+".max");
    /*for (std::vector<set<int> >::iterator a_itr = all_maxima.begin(); a_itr != all_maxima.end(); ++a_itr) {
        for (std::set<int>::iterator itr = a_itr->begin(); itr != a_itr->end(); ++itr) {
            std::cout << *itr << " : " << stabilities[*itr] << " ";
        }
        cout << endl;
    }*/



    cout << "Loading graph into file for bottleneck finding" << endl;
    SmartGraph g;
    graphReader(g, UNIPROJ_HOME+"/data/landscapes/"+basename+"_landscape.lgf").run();

    std::vector<std::map<std::pair<int,int>,float> > all_pairwise_bottlenecks;
    {
        IterableValueMap<SmartGraph, Edge, float> weights(g);
        cout << "Is Connected : " << connected(g) << endl;


        //cout << "Generating Weights" << endl;
        for (int offset = 0; offset < range*num_partitions; offset = offset + num_partitions) {
            float tau = taus[offset/num_partitions];
            cout << "t =  " << tau << endl;
            std::map<std::pair<int,int>,float> pairwise_bottlenecks;
            for (SmartGraph::EdgeIt itr(g); itr != INVALID; ++itr) {
                float u_stab = stabilities[g.id(g.u(itr)) + offset];
                float v_stab = stabilities[g.id(g.v(itr)) + offset];
                float lesser_stab = u_stab < v_stab ? u_stab : v_stab;
                weights.set(itr,lesser_stab);
            }


            set<int> current_set = all_maxima[offset/num_partitions];
            for(set<int>::const_iterator iter1 = current_set.begin(); iter1 != current_set.end(); ++iter1) {
                for(set<int>::const_iterator iter2 = iter1; ++iter2 != current_set.end();)
                {
                    std::cout << *iter1 << " " << *iter2 << " ";
                    SmartGraph::Node source = g.nodeFromId(*iter1);
                    SmartGraph::Node target = g.nodeFromId(*iter2);
                    SmartGraph::EdgeMap<bool> filter(g, true);
                    pairwise_bottlenecks[make_pair(*iter1,*iter2)] = find_bottleneck(g,weights,source,target,filter);
                    cout << pairwise_bottlenecks[make_pair(*iter1,*iter2)] << endl;
                }
            }
            all_pairwise_bottlenecks.push_back(pairwise_bottlenecks);
        }
    }
    save_bottlenecks(all_pairwise_bottlenecks,UNIPROJ_HOME+"/data/bottlenecks/"+basename+".bnp");
    //load_bottlenecks(all_pairwise_bottlenecks,UNIPROJ_HOME+"/data/bottlenecks/"+basename+".bnp");
    {
        std::string filename = UNIPROJ_HOME+"/data/dg/"+basename+".dg";//+energy_string.str();
        ofstream dg_file(filename.c_str());

        for (int offset = 0; offset < range*num_partitions; offset = offset + num_partitions) {
            float tau = taus[offset/num_partitions];
            cout << "t =  " << tau << endl;
            std::map<std::pair<int,int>,float> pairwise_bottlenecks = all_pairwise_bottlenecks[offset/num_partitions];
            set<int> current_set = all_maxima[offset/num_partitions];
            std::map<int,int> super_basins;
            std::set<float> energies;

            for (std::map<std::pair<int,int>,float>::iterator itr = pairwise_bottlenecks.begin(); itr != pairwise_bottlenecks.end(); ++itr) {
                energies.insert(itr->second);
            }

            for (std::set<float>::iterator eitr = energies.begin(); eitr != energies.end(); ++eitr) {
                float energy = *eitr;

            //for (float energy = -1.0; energy < 1.0; energy = energy + 0.1) {
                std::map<int,int> super_basins;
                int i = 0;
                for (std::set<int>::iterator itr = current_set.begin(); itr != current_set.end(); ++itr) {
                    if (stabilities[*itr + offset] > energy) {
                        super_basins[*itr] = i;
                        i++;
                    }
                }
                cout << "energy: " << energy << endl;
                print_map(super_basins);

                if (super_basins.size() > 0) {
                    for(std::map<int,int>::iterator iter1 = super_basins.begin(); iter1 != super_basins.end(); ++iter1) {
                        for(std::map<int,int>::iterator iter2 = iter1; ++iter2 != super_basins.end();) {
                            if (pairwise_bottlenecks[make_pair(iter1->first,iter2->first)] >= energy) {
                                int maxima_1_group = iter1->second;
                                int maxima_2_group = iter2->second;
                                for (std::map<int,int>::iterator sbitr = super_basins.begin(); sbitr != super_basins.end(); ++sbitr) {
                                    if (sbitr->second == maxima_1_group) {
                                        sbitr->second = maxima_2_group;
                                    }
                                }
                            }
                        }
                    }

                    std::ostringstream energy_string;
                    energy_string << energy;
                    reorder_map(super_basins);
                    cout << "saving" << endl;
                    print_map(super_basins);

                    for (std::map<int,int>::iterator itr = super_basins.begin(); itr != super_basins.end(); ++itr) {
                        dg_file << tau << "," << energy << "," << itr->first << "," << itr->second << endl;
                    }
                }
            }
        }
        dg_file.close();
    }

    // TODO TODAY
    // save pairs
    // load pairs and find good resolution for dg graph
    // Implmeent complexity measure for graphs



    //Create new file for each energy level, where line=maxima_id, super_basin_id

    //Save Bottleneck
    //for (std::map<std::pair<int,int>,float>::iterator itr = pairwise_bottlenecks.begin(); itr != pairwise_bottlenecks.end(); ++itr) {

    //}


    // Random Walk
    /*{
        cout << "Loading graph into file for random walk" << endl;
        SmartGraph g;
        graphReader(g, UNIPROJ_HOME+"/data/landscapes/"+basename+"_landscape.lgf").run();
        SmartGraph::NodeMap<int> degrees(g);
        int num_nodes = countNodes(g);

        for (SmartGraph::NodeIt itr(g); itr != INVALID; ++itr) {
            int degree = 0;
            for (SmartGraph::IncEdgeIt e_itr(g,itr); e_itr != INVALID; ++e_itr) {
                degree++;
            }
            degrees[itr] = degree;
        }

        std::vector<float> complexities;

        for (int offset = 0; offset < range*num_partitions; offset = offset + num_partitions) {
            int num_steps = 10000;
            int num_walks = 30;
            float eta_star = 0;
            int num_eta_steps = 200;
            random_walk(30000, 0.0, g, &stabilities, eta_star, degrees, offset);
            //cout << "eta is " << eta_star << endl;

            float max_entropy = 0;
            std::map<float,float> mean_entropies;
            //mean_entropies.insert(mean_entropies.begin(),num_eta_steps+1,0.0);
            for (float eta = 0.0; eta < eta_star; eta = eta + eta_star/num_eta_steps) {
                mean_entropies[eta] = 0.0;
            }

            for (int i = 0; i < num_walks; ++i) {
                for (float eta = 0.0; eta < eta_star; eta = eta + eta_star/num_eta_steps) {
                    float temp_eta_star;
                    std::vector<int> symbols = random_walk(num_steps, eta, g, &stabilities, temp_eta_star, degrees, offset);

                    float current_entropy = calculate_entropy(symbols);
                    mean_entropies[eta] += current_entropy;
                    //cout << current_entropy << endl;

                    //if (current_entropy > max_entropy) {
                    //    max_entropy = current_entropy;
                    //}
                }
            }

            for (float eta = 0.0; eta < eta_star; eta = eta + eta_star/num_eta_steps) {
                mean_entropies[eta] = mean_entropies[eta] / num_walks;
                //cout << mean_entropies[eta] << endl;

                if (mean_entropies[eta] > max_entropy) {
                    max_entropy = mean_entropies[eta];
                }
            }
            //getchar();
            complexities.push_back(max_entropy);
        }
        save_vector(complexities,UNIPROJ_HOME+"/data/complexities/"+basename+".entropy");
    }*/

    //TASKS LEFT
    //CORRELATION
    //DGGRAPH
    //STABILITYGRAPH + MEASURES
    //FIXNEIGHBOURSBUG
    //STABILITY

    // Time series correlation
    // DG Graph
    // Compute All Complexity measures for
        // Stability range for barbell graph
        // Series of random graphs of increasing n or m
        // Varietyy of small graphs for all stability range

    // folder data/test_graphs -- full of .nkes
            // Store DG graph/ DG complexity (line, t,
            // Store auto correlation (line:t, complexity (at t))

    // Hierharchical Clustering of stability
        //TODO From DATA (node:cluster id at given time)
            // Load into memoery and create graph where nodes =  partitions at a certain time
            // and edges between t and t+1 show that particular node is in belongs in boht
            // Could do at every markov time sample Or when partitions change
        // Find measures from this graph
            // Size : V,E
            // between ness, small worldnesesness etc
        // Compar measures between normal patient and schizo

    // Stability analysis at parcellation 90
    // Stability analysis at parcellation 500
    // Plot on smae graph
        // TODO: Convert coords graph into format required for stability code

    delete[] stabilities;
    return 0;
}
