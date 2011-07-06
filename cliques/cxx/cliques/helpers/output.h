#pragma once

#include <iostream>
#include <fstream>

#include <vector>
#include <map>

namespace cliques {

void print_map(std::map<int, int> my_map, bool add_one = false) {
    for (std::map<int, int>::iterator itr = my_map.begin(); itr != my_map.end(); ++itr) {
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
    basins_file.open(filename);

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

template <typename G>
void graph_to_edgelist_file(std::string filename,G &graph) {
    typedef typename G::EdgeIt EdgeIt;
    std::ofstream graph_file;
    graph_file.open(filename);
    for (EdgeIt e(graph); e != lemon::INVALID; ++e) {
        auto n1 = graph.u(e);
        auto n2 = graph.v(e);
        graph_file << graph.id(n1) << " " << graph.id(n2)
                << std::endl;
    }
    graph_file.close();
}

}
