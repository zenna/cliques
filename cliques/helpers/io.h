/* Copyright (c) Zenna Tavares & M. Schaub - zennatavares@gmail.com, 2010-2011 \
 * Input/Output helpers for writing and reading from disk */
#pragma once

#include <iostream>
#include <fstream>
#include <set>

#include <vector>
#include <map>

namespace cliques {

template<typename G>
void read_edgelist(G &graph, std::string filename) {
    std::ifstream maxima_file(filename.c_str());
    std::string line;
    std::string maxima;

    if (!maxima_file.is_open()) {
        std::cout << "couldn't open file" << std::endl;
        exit(1);
    }

    typedef typename G::Node Node;
    std::map<int, Node> id_to_node;
    while (std::getline(maxima_file, line)) {
        std::stringstream lineStream(line);
        std::set<int> current_maxima;

        std::getline(lineStream, maxima, ' ');
        int node1_id = atoi(maxima.c_str());
        std::getline(lineStream, maxima, ' ');
        int node2_id = atoi(maxima.c_str());

        typename std::map<int, Node>::iterator itr = id_to_node.find(node1_id);
        Node node1, node2;

        if (itr == id_to_node.end()) {
            node1 = graph.addNode();
            id_to_node[node1_id] = node1;
        } else {
            node1 = itr->second;
        }
        itr = id_to_node.find(node2_id);
        if (itr == id_to_node.end()) {
            node2 = graph.addNode();
            id_to_node[node2_id] = node2;
        } else {
            node2 = itr->second;
        }
        graph.addEdge(node1, node2);
    }
    maxima_file.close();
}

template<typename G, typename E>
bool read_edgelist_weighted(std::string filename, G &graph, E &weights) {
    // initialise input stream and strings for readout
    std::ifstream maxima_file(filename.c_str());
    std::string line;
    std::string maxima;

    // check if file is open
    if (!maxima_file.is_open()) {
        std::cout << "couldn't open file:" << filename << std::endl;
        return false;
    }

    // define Node class for convenience
    typedef typename G::Node Node;
    // mapping from id to node
    std::map<int, Node> id_to_node;

    //readout contents from maxima_file into string, line by line
    while (std::getline(maxima_file, line)) {

        std::stringstream lineStream(line);
        //readout node id and weights
        std::getline(lineStream, maxima, ' ');
        int node1_id = atoi(maxima.c_str());
        std::getline(lineStream, maxima, ' ');
        int node2_id = atoi(maxima.c_str());
        std::getline(lineStream, maxima, ' ');
        float weight = atof(maxima.c_str());

        typename std::map<int, Node>::iterator itr = id_to_node.find(node1_id);
        Node node1, node2;

        // If the node is not in the map
        // then create node and add to map
        if (itr == id_to_node.end()) {
            node1 = graph.addNode();
            id_to_node[node1_id] = node1;
        } else {
            node1 = itr->second;
        }

        // same for node 2
        itr = id_to_node.find(node2_id);
        if (itr == id_to_node.end()) {
            node2 = graph.addNode();
            id_to_node[node2_id] = node2;
        } else {
            node2 = itr->second;
        }

        //std::cout << "adding " << graph.id(node1) << " - " << graph.id(node2) << std::endl;
        typename G::Edge edge = graph.addEdge(node1, node2);
        weights.set(edge, weight);
    }

    maxima_file.close();
    return true;
}

template<typename G, typename E>
void write_edgelist_weighted(std::string filename, G &graph, E &weights) {
    // initialise input stream and strings for readout
    std::ofstream maxima_file(filename.c_str()); // This was std::ios_base, any reason?
    std::string maxima;

    // check if file is open
    if (!maxima_file.is_open()) {
        std::cout << "couldn't open file" << std::endl;
        exit(1);
    }

    for(typename G::EdgeIt e(graph); e!=lemon::INVALID; ++e){
        int node1 = graph.id(graph.u(e));
        int node2 = graph.id(graph.v(e));
        double weigth = weights[e];

        maxima_file << node1 << " " << node2 << " " << weigth << "\n";
    }

    maxima_file.close();
}

template<typename P>
void read_partitions_file(std::vector<P> &all_partitions, std::string filename) {
    std::ifstream partitions_file(filename.c_str());
    std::string line;
    std::string set;

    if (!partitions_file.is_open()) {
        std::cout << "couldn't open file:" << filename << std::endl;
        exit(1);
    }

    while (std::getline(partitions_file, line)) {
        std::stringstream lineStream(line);
        std::set<int> current_maxima;

        std::vector<int> raw_partition;
        while(std::getline(lineStream, set, ' ')) {
            raw_partition.push_back(atoi(set.c_str()));
        }
        P nth_partition(raw_partition);
        all_partitions.push_back(nth_partition);
    }
    partitions_file.close();
}

template <typename T>
void print_map(typename std::map<int, T> my_map, bool add_one = false) {
    for (typename std::map<int, T>::iterator itr = my_map.begin(); itr != my_map.end(); ++itr) {
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

template <typename T>
void print_2d_vector(std::vector<std::vector<T> > my_vector) {
    typename std::vector<std::vector<T> >::iterator itr;
    for (itr = my_vector.begin(); itr != my_vector.end(); ++itr) {
        typename std::vector<T>::iterator new_itr;
        for (new_itr = itr->begin(); new_itr != itr->end(); ++new_itr) {
            std::cout << (*new_itr) << " ";
        }
        std::cout << "\n";
    }
}

template<class T>
void print_collection(T collection) {
    for (class T::iterator itr = collection.begin(); itr != collection.end(); ++itr) {
        std::cout << *itr << ", ";
    }
    std::cout << std::endl;
}

template<class T>
void print_collection(T collection, int new_line) {
    int i =0;
    for (class T::iterator itr = collection.begin(); itr != collection.end(); ++itr) {
        if (i % new_line == 0) {
            std::cout << "\n";
        }
        std::cout << *itr << ", ";
        ++i;

    }
    std::cout << std::endl;
}


//Since variadic templates are recursive, must have a base case
void output() { std::cout << '\n'; }

//Output function to output any type without type specifiers like printf() family
template <typename T, typename ...P>
void output(T t, P ...p)
{
  std::cout << t << ' ';
  if (sizeof...(p)) { output(p...); }
  else { std::cout << '\n'; }
}

template<typename P>
void print_partition(P &partition) {
    for (typename P::PartIterator pitr = partition.begin(); pitr
            != partition.end(); ++pitr) {
        std::cout << "[";
        for (typename P::NodeIterator nitr = pitr.begin(); nitr != pitr.end(); ++nitr) {
            std::cout << " " << *nitr << " ";
        }
        std::cout << "]" << std::endl;
    }
}

template<typename P>
void print_partition_list(P &partition) {
    int length = partition.element_count();
    for (int i = 0; i < length; i++) {
        std::cout << i << "->" << partition.find_set(i) << std::endl;
    }
}

template<typename P>
void print_partition_line(P &partition) {
    int length = partition.element_count();
    for (int i = 0; i < length; i++) {
        std::cout << partition.find_set(i) << " ";
    }
    std::cout << std::endl;
}

void basins_to_file(std::string filename,
        std::vector<std::map<int, std::map<int, double> >> all_basins,
        std::vector<double> markov_times) {
    std::ofstream basins_file;
    basins_file.open(filename.c_str());

    for (unsigned int i = 0; i < markov_times.size(); ++i) {
        double time = markov_times[i];
        auto basins = all_basins[i];
        for (auto itr = basins.begin(); itr != basins.end(); ++itr) {
            int basin_id = itr->first;
            basins_file << time << " " << basin_id << " ";
            cliques::output(time, itr->second.size());

            for (auto b = itr->second.begin(); b != itr->second.end(); ++b) {
                basins_file << b->first << " " << b->second << " ";
            }
            basins_file << std::endl;
        }
    }
    basins_file.close();
}

/**
 @brief  Write partitions from a container into a file. (Template)

 This functions iterates over a given container and writes the partitions into 
 a file. Each line stands for one partition, each column stands for one node 
 with its corresponding community Id.

 @param[in]  filename          file to be written
 @param[in]  all_partitions    container with partitions

 */
template <typename P>
void partitions_to_file(std::string filename,
        P & all_partitions) {
    // init streams
    std::ofstream partitions_file;
    partitions_file.open(filename);
    
    // iterate over container and write partitions in file
    for (auto itr = all_partitions.begin();
        itr != all_partitions.end();++itr) {

        int length = itr->element_count();
        for (int i = 0; i < length; i++) {
            partitions_file << itr->find_set(i) << " ";
        }

        partitions_file << std::endl;
    }
    partitions_file.close();
}


template <typename G>
void graph_to_edgelist_file(std::string filename,G &graph) {
    typedef typename G::EdgeIt EdgeIt;
    std::ofstream graph_file;
    graph_file.open(filename.c_str());
    for (EdgeIt e(graph); e != lemon::INVALID; ++e) {
        auto n1 = graph.u(e);
        auto n2 = graph.v(e);
        graph_file << graph.id(n1) << " " << graph.id(n2)
                << std::endl;
    }
    graph_file.close();
}

}
