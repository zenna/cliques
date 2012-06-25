/* Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010-2011 */
/* Algorithms for computing basins and basin properties */
#pragma once

namespace clq {

/**
 @brief  Compute the probabalistic basins defined by kernginal lin
 algorithm
 */
template<typename G, typename M, typename QF, typename QFDIFF, typename MAP,
        typename COMM>
std::map<int, std::map<int, double> > compute_kernighan_lin_basins(G &graph,
        M &weights, QF compute_quality, QFDIFF compute_quality_diff, MAP &map,
        G &space, double time, COMM &all_partitions) {

    std::map<int, std::map<int, int> > node_to_basin_to_count;
    typedef clq::VectorPartition VecPartition;
    const int TOTAL_COUNT_INDEX = -1; // Magic number to store total count
    int num_nodes = lemon::countNodes(graph);

    for (auto start_partition = all_partitions.begin(); start_partition
            != all_partitions.end(); ++start_partition) {

        std::vector<VecPartition> path_partitions;
        VecPartition output_partition(num_nodes);

        clq::refine_partition_kernighan_lin_hijack(graph, weights,
                compute_quality, compute_quality_diff, path_partitions,
                *start_partition, output_partition);

        output_partition.normalise_ids();

        //Convert path_partitions into path_nodes
        std::vector<int> path_nodes;
        //        clq::output("paths");
        for (auto partition = path_partitions.begin(); partition
                != path_partitions.end(); ++partition) {
            partition->normalise_ids();
            //            clq::print_partition_line(*partition);
            auto p = map.left.find(*partition);
            if (p != map.left.end()) {
                auto node_id = space.id(p->second);
                path_nodes.push_back(node_id);
            }
        }

        //        clq::output("start then output:");
        //        clq::print_partition_line(*start_partition);
        //        clq::print_partition_line(output_partition);

        auto p = map.left.find(output_partition);
        if (p != map.left.end()) {
            int optimum = space.id(p->second);
            path_nodes.push_back(space.id(map.left.at(*start_partition)));

            for (auto node = path_nodes.begin(); node != path_nodes.end(); ++node) {
                node_to_basin_to_count[*node][optimum]++;
                node_to_basin_to_count[*node][TOTAL_COUNT_INDEX]++;

                //                clq::output(*node,optimum,node_to_basin_to_count[*node][optimum],node_to_basin_to_count[*node][TOTAL_COUNT_INDEX]);
            }
        } else {
            clq::output("maxima is disconnected partition");
        }

    }

    std::map<int, std::map<int, double> > basin;

    // Convert into basin format
    for (auto node = node_to_basin_to_count.begin(); node
            != node_to_basin_to_count.end(); ++node) {
        int node_id = node->first;
        int total_count = node->second[TOTAL_COUNT_INDEX];

        for (auto basin_itr = node->second.begin(); basin_itr
                != node->second.end(); ++basin_itr) {
            int basin_id = basin_itr->first;
            if (basin_id != -1) {
                int basin_count = basin_itr->second;
                basin[basin_id][node_id] = basin_count / float(total_count);
            }
        }
    }

    return basin;
}

/**
 @brief  Convert a set of partitions into a graph
 */
template<typename G, typename DG, typename M>
void graph_to_transition_digraph(G &graph, std::vector<double> qualities,
        DG &transition_graph, M &transition_weights) {
    typedef typename G::NodeIt NodeIt;
    typedef typename G::IncEdgeIt IncEdgeIt;
    typedef typename G::Node Node;
    typedef typename DG::Arc Arc;

    int num_nodes = lemon::countNodes(graph);
    for (int i = 0; i < num_nodes; ++i) {
        transition_graph.addNode();
    }
    for (NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
        int base_id = graph.id(n1);

        double current_quality = qualities[base_id];
        double total_weight = 0.0;
        std::vector<Arc> unormalised_arcs;
        for (IncEdgeIt e(graph, n1); e != lemon::INVALID; ++e) {
            Node running = graph.oppositeNode(n1, e);
            int running_id = graph.id(running);
            double running_quality = qualities[running_id];
            if (running_quality > current_quality) {
                double weight_difference = running_quality - current_quality;
                total_weight += weight_difference;
                Arc new_arc = transition_graph.addArc(
                        transition_graph.nodeFromId(base_id),
                        transition_graph.nodeFromId(running_id));
                unormalised_arcs.push_back(new_arc);
                transition_weights.set(new_arc, weight_difference);
            }
        }
        for (typename std::vector<Arc>::iterator arc = unormalised_arcs.begin(); arc
                != unormalised_arcs.end(); ++arc) {
            double normalised_transition_weight = transition_weights[*arc]
                    / total_weight;
            transition_weights.set(*arc, normalised_transition_weight);
        }
    }
}

/**
 @brief  Given a container of (generic) configurations

 @param[in]  all_configurations iterable set of sampled configurations
 @param[out]  basins            Optima -> (NodeGoingToOptima -> Probability of going there)
 */
template<typename C, typename cC>
std::map<int, std::map<int, double> >
sample_probabalistic_basins( 
        cC &all_configs,
        std::function<C (C)> optimisation_func,
        const int num_samples,
        std::function<int (C)> get_id) {

    std::map<int, std::map<int, int> > config_to_basin_to_count;
    cC basins;

    // typedef typename clq::VectorPartition VecPart;
    // typedef typename std::unordered_set<VecPart, clq::partition_hash, clq::partition_equal> partition_set;

    std::unordered_map<C, std::unordered_map<C, int, clq::partition_hash, clq::partition_equal>,
        clq::partition_hash, clq::partition_equal > config_to_basin_to_count_x;
    clq::output("size of all configs before ", all_configs.size());

    // For all configs, optimise to optimum num_samples times
    int jj = 0;
    for (C const &config : all_configs) {
        clq::output("Sampling from i=",jj, "from total of", all_configs.size());
        int config_id = jj;//get_id(config);
        for (int i = 0;i<num_samples;++i) {
            // clq::output("sample ",i," of ",num_samples);
            // clq::output("input partition is");
            // clq::print_partition_line(config);
            C basin = optimisation_func(config); // INVALID READ HERE
            // clq::output("output is");
            // clq::print_partition_line(basin);

            // If the basin is not found then add it
            // FIXME: make all_configs const, and return a new container
            // Make independent of vector type (remove push_back)
            config_to_basin_to_count_x[config][basin]++;

            int basin_id = get_id(basin);
            if (basin_id == -1 )  { // magic number for not found
                basins.push_back(basin); // AND HERE
                clq::output("basin not in original list, adding");
            }

            // config_to_basin_to_count[config_id][get_id(basin)]++;
        }
        ++jj;
    };

    // Then add basins to list 
    for (C const &basin : basins) {
        all_configs.push_back(basin);
    }

    clq::output("size of all configs after ", all_configs.size());

    std::map<int, std::map<int, double> > basin_to_config_to_prob;    

    // For very nodes n, normalise counts to find prob of going to each basin
    for (auto &config_to_basin_to_count1 : config_to_basin_to_count_x) {
        int config_id = get_id(config_to_basin_to_count1.first);
        auto basin_to_count = config_to_basin_to_count1.second;

        for (auto basin : basin_to_count) {
            int basin_id = get_id(basin.first);
            int count = basin.second;
            basin_to_config_to_prob[basin_id][config_id] = (double)count / (double)num_samples;
        }
    }

    return basin_to_config_to_prob;
}

/**
 @brief  Compute probabilistic basins using message passing algorithm

 */
template<typename G>
std::map<int, std::map<int, double> > compute_mp_probabalistic_basins(
        G &graph, std::vector<double> qualities) {

    typedef typename lemon::SmartDigraph::NodeIt NodeIt;
    typedef typename lemon::SmartDigraph::Node Node;
    typedef typename lemon::SmartDigraph::OutArcIt OutArcIt;

    lemon::SmartDigraph transition_graph;
    lemon::SmartDigraph::ArcMap<double> transition_weights(transition_graph);
    graph_to_transition_digraph(graph, qualities, transition_graph,
            transition_weights);

    std::map<int, std::map<int, double> > all_basins;
    std::set<int> maxima;
    std::set<Node> maxima_nodes;
    lemon::SmartDigraph::NodeMap<int> real_outlinks_queue(transition_graph, 0);

    //find all maxima and insert into set
    for (NodeIt node(transition_graph); node != lemon::INVALID; ++node) {
        int num_better_nodes = 0;
        for (OutArcIt arc(transition_graph, node); arc != lemon::INVALID; ++arc) {
            num_better_nodes++;
        }
        real_outlinks_queue[node] = num_better_nodes;
        // I.e. only if this node is a maxima
        if (num_better_nodes == 0) {
            clq::output("maxima");
            int node_id = transition_graph.id(node);
            maxima.insert(node_id);
            maxima_nodes.insert(node);
        }

    }

    for (std::set<int>::iterator itr = maxima.begin(); itr != maxima.end(); ++itr) {
        // for each maximum do message passing
        lemon::SmartDigraph::NodeMap<double> node_values(transition_graph, 0);
        int max_node_id = *itr;
        Node max_node = transition_graph.nodeFromId(max_node_id);
        node_values[max_node] = 1;
        std::set<Node> candidates(maxima_nodes);
        lemon::SmartDigraph::ArcMap<double> weights_to_max(transition_graph);
        lemon::mapCopy(transition_graph, transition_weights, weights_to_max);

        lemon::SmartDigraph::NodeMap<int> outlinks_queue(transition_graph);
        lemon::mapCopy(transition_graph, real_outlinks_queue, outlinks_queue);

        // as long as there are "unconverged" Nodes
        while (!candidates.empty()) {

            //loop over Nodes to be considered
            for (std::set<Node>::iterator candidate = candidates.begin(); candidate
                    != candidates.end(); ++candidate) {
                Node candidate_node = *candidate;

                // sum up outgoing link weight..
                for (lemon::SmartDigraph::OutArcIt outgoing_link(
                        transition_graph, candidate_node); outgoing_link
                        != lemon::INVALID; ++outgoing_link) {
                    node_values[candidate_node]
                            += weights_to_max[outgoing_link];
                }

                // multiply node value "down"..
                for (lemon::SmartDigraph::InArcIt incoming_link(
                        transition_graph, candidate_node); incoming_link
                        != lemon::INVALID; ++incoming_link) {

                    weights_to_max[incoming_link]
                            *= node_values[candidate_node];
                    Node downstream_node = transition_graph.source(
                            incoming_link);
                    outlinks_queue[downstream_node] -= 1;

                    // add potential nodes to candidate list
                    if (outlinks_queue[downstream_node] == 0) {
                        candidates.insert(downstream_node);
                    }

                }
                // remove node that passed on message
                candidates.erase(candidate_node);
            }
        }

        std::map<int, double> basin_to_probabilities;
        for (NodeIt node(transition_graph); node != lemon::INVALID; ++node) {
            int node_id = transition_graph.id(node);
            if (node_values[node]) {
                basin_to_probabilities[node_id] = node_values[node];
            }
        }

        all_basins[max_node_id] = basin_to_probabilities;
    }
    return all_basins;
}

}