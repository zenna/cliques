#include <vector>
#include <iostream>

#include <algorithm>
#include <vector>
#include <map>

#include <lemon/smart_graph.h>
#include <boost/program_options.hpp>

#include <cliques/helpers.h>
#include <cliques/algorithms/all_partitions.h>
#include <cliques/algorithms/all_communities.h>
#include <cliques/algorithms/stability.h>
#include <cliques/algorithms/louvain.h>
#include <cliques/algorithms/space.h>
#include <cliques/algorithms/maxima.h>
#include <cliques/algorithms/kernighan_lin.h>
#include <cliques/algorithms/internals/internals.h>
#include <cliques/nldr/nldr.h>
#include <cliques/structures/make_graphs.h>
#include <cliques/structures/vector_partition.h>
#include <cliques/algorithms/aglob.h>


namespace po = boost::program_options;

template<typename G, typename M>
void parse_arguments(int ac, char *av[], G &graph, M &weights,
        int &num_samples, int &num_dim, std::string &filename_prefix) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")
    		("graph,G",po::value<std::string>(), "input graph")
    		("num-samples,S",po::value<int>(), "number of samples")
    		("dimensions,d",po::value<int>(), "number of dimensions")
            ("prefix,x",po::value<std::string>(),"filename prefix");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        //return 1; - exit
    }

    if (vm.count("num-samples")) {
        num_samples = vm["num-samples"].as<int> ();
    }

    if (vm.count("dimensions")) {
        num_dim = vm["dimensions"].as<int> ();
    }

    if (vm.count("prefix")) {
        filename_prefix = vm["prefix"].as<std::string> ();
    }
    else {
    	filename_prefix = "out";
    }

    cliques::output("making graph");
    if (vm.count("graph")) {
        std::string filename = vm["graph"].as<std::string> ();
        cliques::read_edgelist_weighted(filename, graph, weights);
    } else {
        //cliques::make_path_graph(graph, 7, weights);
        //      cliques::make_ring_graph(graph, 12, weights);
        cliques::make_complete_graph(graph, 7, weights);
    }
}

int main(int ac, char* av[]) {
    typedef cliques::VectorPartition VecPartition;
    typedef boost::unordered_set<VecPartition, cliques::partition_hash,
            cliques::partition_equal> VecPartitionSet;

    lemon::SmartGraph orange_graph;
    std::string filename_prefix;
    lemon::SmartGraph::EdgeMap<double> weights(orange_graph);
    int num_samples = 100000;
    int num_dim = 3;
    parse_arguments(ac, av, orange_graph, weights, num_samples, num_dim, filename_prefix);
    bool find_partitions, create_space, find_stabs, find_basins, find_dists;
    find_partitions = create_space = find_stabs = find_basins = find_dists = false;

    cliques::output("Find Connected Communities");
    auto a = cliques::find_connected_communities(orange_graph);
//    cliques::print_2d_vector(a);
//    cliques::output(a.size());
//    cliques::output("neighbours of ");
//    cliques::print_collection(a[0]);
//    std::vector<double> markov_times = {1.0};
//    auto b = cliques::find_community_neighbours(orange_graph, a[0]);
//    cliques::find_weighted_linearised_stability func(markov_times);
//
//    cliques::print_2d_vector(b);
//    for (auto comm = a.begin(); comm != a.end(); ++comm) {
//        double stability;
////        cliques::print_collection(*comm);
//        bool is_max = cliques::is_community_maxima(orange_graph,  weights, *comm, func, stability);
//        if (is_max) {
//            cliques::output("maxima:", stability);
//            cliques::print_collection(*comm);
//        }
//    }

    cliques::output("Finding Connected Partitions");
    cliques::NoLogging no_logging;
    cliques::Logging<VecPartition> log_all;
    VecPartitionSet all_partitions;
    cliques::find_connected_partitions(orange_graph, all_partitions, no_logging);
    cliques::output("complete size:", all_partitions.size());
//
//    cliques::output("Creating space graph");
//    lemon::SmartGraph space;
//    auto map = cliques::create_space(orange_graph, all_partitions, space);
//    lemon::SmartGraph::EdgeMap<float> space_weights(space);
//    cliques::make_weights_from_edges(space, space_weights);
//
//    cliques::output("Finding stabilities");
//    std::ofstream stabs_file;
//    stabs_file.open(filename_prefix + "_energy.mat");
//    std::vector<double> markov_times;
//    for (double t = 0.01; t < 5.0; t = t * 1.005) {
//        markov_times.push_back(t);
//    }
//    cliques::output(markov_times.size());
//    cliques::find_weighted_linearised_stability func(markov_times);
//    for (unsigned int i = 0; i < markov_times.size(); ++i) {
//        stabs_file << markov_times[i] << " ";
//        for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
//            cliques::LinearisedInternals internals(orange_graph, weights, *itr);
//            double stability = func(internals, i);
//            stabs_file << stability << " ";
//        }
//        if (i + 1 != markov_times.size()) {
//            stabs_file << std::endl;
//        }
//    }
//    stabs_file.close();
//
//    cliques::output("Finding Maxima");
//    VecPartitionSet all_maxima;
//    double current_markov_time = markov_times[0];
//    cliques::linearised_stability_gain_louvain diff_func(current_markov_time);
//    cliques::sample_maxima(orange_graph, weights, func, diff_func, all_maxima,
//            all_partitions, no_logging);
//
//    std::ofstream graph_file;
//    graph_file.open(filename_prefix + "_graph_edgelist.edj");
//    for (lemon::SmartGraph::EdgeIt e(orange_graph); e != lemon::INVALID; ++e) {
//        auto n1 = orange_graph.u(e);
//        auto n2 = orange_graph.v(e);
//        graph_file << orange_graph.id(n1) << " " << orange_graph.id(n2)
//                << std::endl;
//    }
//    graph_file.close();
//
////    cliques::output("Kernighan Lin Basin of attractions");
////    std::ofstream basin_file;
////    basin_file.open("basin.mat");
////    // Format: time optima_node basin_node basin_value basin_node basin_value
////    int num_nodes = lemon::countNodes(orange_graph);
////    for (unsigned int i = 0; i < markov_times.size(); ++i) {
////        std::map<int, std::vector<int> > optima_to_basin;
////        for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
////            cliques::VectorPartition optima(num_nodes);
////            std::vector<double> mkov_times;
////            mkov_times.push_back(markov_times[i]);
////            double mkov = markov_times[i];
////            cliques::refine_partition_kernighan_lin(orange_graph, weights,
////                    cliques::find_weighted_linearised_stability(mkov_times),
////                    cliques::linearised_stability_gain_louvain(mkov), *itr,
////                    optima);
////            optima.normalise_ids();
////            int optima_id =  orange_graph.id(map.left.at(optima));
////            int partition_id = orange_graph.id(map.left.at(*itr));
////            optima_to_basin[optima_id].push_back(partition_id);
////        }
////        cliques::output("num_basins", optima_to_basin.size());
////        for (auto itr = optima_to_basin.begin(); itr != optima_to_basin.end(); ++itr) {
////            lemon::SmartGraph::Node n = orange_graph.nodeFromId(itr->first);
////            cliques::print_partition_line(map.right.at(n));
////            cliques::output("time: ", markov_times[i], "num in basin", itr->second.size());
////            basin_file << markov_times[i] << " " << itr->first << " ";
////            for (auto b_itr = itr->second.begin(); b_itr != itr->second.end(); ++b_itr) {
////                basin_file << *b_itr << " " << "1.0 ";
////            }
////            basin_file << std::endl;
////        }
////    }
////    basin_file.close();
//
////    cliques::output("Louvain Basin of attractions");
////    std::ofstream basin_file;
////    basin_file.open("basin.mat");
////    // Format: time optima_node basin_node basin_value basin_node basin_value
////    int num_nodes = lemon::countNodes(orange_graph);
////    for (unsigned int i = 0; i < markov_times.size(); ++i) {
////        std::map<int, std::vector<int> > optima_to_basin;
////        for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
////            cliques::VectorPartition optima(num_nodes);
////            std::vector<double> mkov_times;
////            mkov_times.push_back(markov_times[i]);
////            double mkov = markov_times[i];
////            std::vector<VecPartition> optimal_partitions;
////            cliques::find_optimal_partition_louvain_with_gain(orange_graph,
////                    weights, cliques::find_weighted_linearised_stability(mkov_times),
////                    cliques::linearised_stability_gain_louvain(mkov),
////                    *itr,
////                    optimal_partitions, no_logging);
////
////            optima = optimal_partitions.back();
////
////            optima.normalise_ids();
////            cliques::print_partition_line(*itr);
////            cliques::print_partition_line(optima);
////            int optima_id =  orange_graph.id(map.left.at(optima));
////            cliques::output("between");
////            int partition_id = orange_graph.id(map.left.at(*itr));
////            optima_to_basin[optima_id].push_back(partition_id);
////        }
////        cliques::output("num_basins", optima_to_basin.size());
////        for (auto itr = optima_to_basin.begin(); itr != optima_to_basin.end(); ++itr) {
////            lemon::SmartGraph::Node n = orange_graph.nodeFromId(itr->first);
//////            cliques::print_partition_line(map.right.at(n));
//////            cliques::output("time: ", markov_times[i], "num in basin", itr->second.size());
////            basin_file << markov_times[i] << " " << itr->first << " ";
////            for (auto b_itr = itr->second.begin(); b_itr != itr->second.end(); ++b_itr) {
////                basin_file << *b_itr << " " << "1.0 ";
////            }
////            basin_file << std::endl;
////        }
////    }
////    basin_file.close();
//
//    cliques::output("S Basin of attractions");
//    std::ofstream basin_file;
//    basin_file.open(filename_prefix + "_greedy_basins.mat");
//    // Format: time optima_node basin_node basin_value basin_node basin_value
//    int num_nodes = lemon::countNodes(orange_graph);
//    for (unsigned int i = 0; i < markov_times.size(); ++i) {
//        std::map<int, std::vector<int> > optima_to_basin;
//        for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
//            cliques::VectorPartition optima(num_nodes);
//            std::vector<double> mkov_times;
//            mkov_times.push_back(markov_times[i]);
//            double mkov = markov_times[i];
//
//            VecPartitionSet temp_set;
//            temp_set.insert(*itr);
//
//            VecPartitionSet all_maxima;
//            double current_markov_time = markov_times[0];
//            cliques::linearised_stability_gain_louvain diff_func(current_markov_time);
//            cliques::sample_maxima(orange_graph, weights, cliques::find_weighted_linearised_stability(mkov_times),
//                    cliques::linearised_stability_gain_louvain(mkov), all_maxima,
//                    temp_set, no_logging);
//
//            optima = *(all_maxima.begin());
//
//            optima.normalise_ids();
//            int optima_id =  orange_graph.id(map.left.at(optima));
//            int partition_id = orange_graph.id(map.left.at(*itr));
//            optima_to_basin[optima_id].push_back(partition_id);
//        }
//        cliques::output("num_basins", optima_to_basin.size());
//        for (auto itr = optima_to_basin.begin(); itr != optima_to_basin.end(); ++itr) {
//            lemon::SmartGraph::Node n = orange_graph.nodeFromId(itr->first);
////            cliques::print_partition_line(map.right.at(n));
////            cliques::output("time: ", markov_times[i], "num in basin", itr->second.size());
//            basin_file << markov_times[i] << " " << itr->first << " ";
//            for (auto b_itr = itr->second.begin(); b_itr != itr->second.end(); ++b_itr) {
//                basin_file << *b_itr << " " << "1.0 ";
//            }
//            basin_file << std::endl;
//        }
//    }
//    basin_file.close();

//    for (auto p = all_maxima.begin(); p != all_maxima.end(); ++p) {
//        cliques::print_partition_line(*p);
//    }

//
//    std::ofstream space_file;
//    space_file.open(filename_prefix + "_space_edgelist.edj");
//    for (lemon::SmartGraph::EdgeIt e(space); e != lemon::INVALID; ++e) {
//        auto n1 = space.u(e);
//        auto n2 = space.v(e);
//        space_file << space.id(n1) << " " << space.id(n2)
//                << std::endl;
//    }
//    space_file.close();

//    cliques::output("Finding distances");
//    //auto X = cliques::find_geodesic_dists(space, landmark_nodes, space_weights);
//    auto X = cliques::find_edit_dists(all_partitions);
//
//    // From edit matrix: find only ones. output into two_d matrix
//    cliques::output(X.n_cols, X.n_rows);
//    std::vector<std::vector<int> > edges;
//    for (unsigned int i = 0; i < X.n_rows; ++i) {
//        for (unsigned int j = i + 1; j < X.n_cols; ++j) {
//            if (X(i, j) == 1) {
//                std::vector<int> edge;
//                edge.push_back(i);
//                edge.push_back(j);
//                edges.push_back(edge);
//            }
//        }
//    }
//    arma::umat edges_mat(edges.size(), 2);
//    int i = 0;
//    for (auto itr = edges.begin(); itr != edges.end(); ++itr) {
//        edges_mat(i, 0) = (*itr)[0];
//        edges_mat(i, 1) = (*itr)[1];
//        ++i;
//    }
//    edges_mat.save(filename_prefix + "_landscape_edgelist.edj", arma::raw_ascii);
////
////    cliques::output("finding embedding");
////    auto L = cliques::embed_mds(X, num_dim);
////    arma::mat L_t = arma::trans(L);
////    L_t.save(filename_prefix + "_coords.mat", arma::raw_ascii);
////
////    auto D_y = cliques::euclid_pairwise_dists(L_t);
////    cliques::output("residual variance", cliques::residual_variance(X, D_y));
//
    std::ofstream vector_file;
    vector_file.open(filename_prefix + "_partitions.mat");
    for (auto itr = all_partitions.begin(); itr != all_partitions.end(); ++itr) {
        int length = itr->element_count();
        for (int i = 0; i < length; i++) {
            vector_file << itr->find_set(i) << " ";
        }
        vector_file << std::endl;
    }

    //    cliques::output("number of nodes", lemon::countNodes(space));
    //    cliques::output("number of edges", lemon::countEdges(space));
    return 0;
}
