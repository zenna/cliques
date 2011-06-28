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
    lemon::SmartGraph::EdgeMap<float> weights(orange_graph);
    int num_samples = 100000;
    int num_dim = 3;
    parse_arguments(ac, av, orange_graph, weights, num_samples, num_dim, filename_prefix);


    //[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    //view.js:1080 0 44 44
    //view.js:114[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0]
//
//    [1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0]
//    view.js:1080 2 32 544
//    view.js:114[0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0]
//

//    std::vector<int> c1 = {1,};
//    std::vector<int> c2 = {0,1,2,3,4,6,7,8,9};
//    std::vector<int> buffer(12);
//    double dist = cliques::find_community_dist(orange_graph, weights, c1, c2, buffer);
//    cliques::output("test_dist:", cliques::find_community_dist(orange_graph, weights, c1, c2, buffer));

    cliques::output("Find Connected Communities");
    auto communities = cliques::find_connected_communities(orange_graph);
    cliques::print_2d_vector(communities);
//    cliques::output(communities.size(), "connected communities");
//    cliques::output("neighbours of ");
//    cliques::print_collection(a[0]);
    std::vector<double> markov_times;
    for (double t = 0.01; t < 5.0; t = t * 1.005) {
        markov_times.push_back(t);
    }
//    auto communities = cliques::find_community_neighbours(orange_graph, a[0]);
    cliques::find_weighted_linearised_stability func(markov_times);

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

//    cliques::output("Creating space graph");
//    lemon::SmartGraph space;
//    auto map = cliques::create_space(orange_graph, all_partitions, space);
//    lemon::SmartGraph::EdgeMap<float> space_weights(space);
//    cliques::make_weights_from_edges(space, space_weights);
////

    cliques::output("Finding stabilities");
    std::ofstream stabs_file;
    stabs_file.open(filename_prefix + "_energy.mat");
    cliques::output(markov_times.size());

    for (unsigned int i = 0; i < markov_times.size(); ++i) {
        stabs_file << markov_times[i] << " ";
        for (auto itr = communities.begin(); itr != communities.end(); ++itr) {
        	cliques::VectorPartition p = cliques::community_to_partition(orange_graph, *itr, 0);
            cliques::LinearisedInternals internals(orange_graph, weights, p);
            double stability = func(internals, 1, i);
            stabs_file << stability << " ";
        }
        if (i + 1 != markov_times.size()) {
            stabs_file << std::endl;
        }
    }
    stabs_file.close();

    std::ofstream graph_file;
    graph_file.open(filename_prefix + "_graph_edgelist.edj");
    for (lemon::SmartGraph::EdgeIt e(orange_graph); e != lemon::INVALID; ++e) {
        auto n1 = orange_graph.u(e);
        auto n2 = orange_graph.v(e);
        graph_file << orange_graph.id(n1) << " " << orange_graph.id(n2)
                << std::endl;
    }
    graph_file.close();

    cliques::output("Finding distances");
    auto X = cliques::find_community_edit_dists(orange_graph,communities);

    // From edit matrix: find only ones. output into two_d matrix
    cliques::output(X.n_cols, X.n_rows);
    std::vector<std::vector<int> > edges;
    for (unsigned int i = 0; i < X.n_rows; ++i) {
        for (unsigned int j = i + 1; j < X.n_cols; ++j) {
            if (X(i, j) == 1) {
                std::vector<int> edge;
                edge.push_back(i);
                edge.push_back(j);
                edges.push_back(edge);
            }
        }
    }
    arma::umat edges_mat(edges.size(), 2);
    int i = 0;
    for (auto itr = edges.begin(); itr != edges.end(); ++itr) {
        edges_mat(i, 0) = (*itr)[0];
        edges_mat(i, 1) = (*itr)[1];
        ++i;
    }
    edges_mat.save(filename_prefix + "_landscape_edgelist.edj", arma::raw_ascii);

    cliques::output("finding embedding");
    auto L = cliques::embed_mds(X, num_dim);
    arma::mat L_t = arma::trans(L);
    L_t.save(filename_prefix + "_coords.mat", arma::raw_ascii);

    auto D_y = cliques::euclid_pairwise_dists(L_t);
    cliques::output("residual variance", cliques::residual_variance(X, D_y));

    std::ofstream vector_file;
    vector_file.open(filename_prefix + "_partitions.mat");
    for (auto itr = communities.begin(); itr != communities.end(); ++itr) {
    	cliques::VectorPartition p = cliques::community_to_partition(orange_graph,*itr,0);
        int length = p.element_count();
        for (int i = 0; i < length; i++) {
            vector_file << p.find_set(i) << " ";
        }
        vector_file << std::endl;
    }

    //    cliques::output("number of nodes", lemon::countNodes(space));
    //    cliques::output("number of edges", lemon::countEdges(space));
    return 0;
}
