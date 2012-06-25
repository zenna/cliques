//  Copyright (c) Zenna Tavares - zennatavares@gmail.com, 2010- 2011 
// /* Algorithms for computing basins and basin properties */
// #pragma once

// #include <limits>
// #include <vector>
// #include <set>

// #include <lemon/concepts/graph.h>

// #include <cliques/quality_functions/internals/internals.h>

// namespace clq {
// /**
//  @brief  Finds the maxima of a landscape graph through steepest ascent

//  This is an exhaustive search to find maximal points on a landscape
//  when a landscape is represented as a graph

//  @param[in] graph the input graph
//  @param[in] stabilities pointer (array) of node 'energies'
//  @tparam G graph type
//  */
// template<typename G>
// std::set<int> find_maxima(G &graph, double *stabilities) {
//     std::set<int> maxima;
//     int num_iterations = 0;
//     for (typename G::NodeIt n(graph); n != lemon::INVALID; ++n) {
//         // Check part_itr hasn't already been visited
//         if (num_iterations % 100 == 0) {
//             //std::cout << "the number of iterations is" << num_iterations << std::endl;
//             //std::cout << "so far number of maxima is " << maxima.size() << std::endl;
//             //std::cout << "size of visited is " << visited.size() << std::endl;
//             //std::cout << "num skipped " << num_skipped << std::endl;
//         }
//         typename G::Node best_neighbour = n;
//         while (1) {
//             bool has_improved = false;
//             double best_score = stabilities[graph.id(best_neighbour)];
//             for (typename G::OutArcIt a(graph, best_neighbour); a
//                     != lemon::INVALID; ++a) {
//                 if (best_score < stabilities[graph.id(graph.target(a))]) {
//                     best_score = stabilities[graph.id(graph.target(a))];
//                     best_neighbour = graph.target(a);
//                     has_improved = true;
//                 }
//             }

//             if (has_improved == false) {
//                 maxima.insert(graph.id(best_neighbour));
//                 break;
//             }
//         }
//         num_iterations++;
//     }
//     return maxima;
// }

// /**
//  @brief  Samples the maxima of a landscape graph through steepest ascent

//  This is an sampled search to find maximal points on a landscape
//  */
// template<typename T, typename W, typename QF, typename sP, typename Logger>
// void sample_maxima(T &graph, W &weights, QF &compute_quality, double time, sP &maxima,
//         sP &sampled_partitions, Logger log) {

//     typedef typename sP::value_type Partition_Type;
//     typedef sP partition_set;

//     Partition_Type partition(lemon::countNodes(graph));
//     partition.initialise_as_singletons();

//     for (auto set_itr = sampled_partitions.begin(); set_itr
//             != sampled_partitions.end(); ++set_itr) {
//         auto best_neighbour = *set_itr;
//         while (true) {
//             bool has_improved = false;
//             partition_set neighs = clq::find_neighbours(graph, best_neighbour);
//             double best_quality = compute_quality(best_neighbour, time);

//             for (auto neigh_itr = neighs.begin(); neigh_itr != neighs.end(); ++neigh_itr) {
//                 double neigh_quality = compute_quality(*neigh_itr, time);
//                 if (best_quality < neigh_quality) {
//                     best_quality = neigh_quality;
//                     best_neighbour = *neigh_itr;
//                     has_improved = true;
//                 }
//             }
//             if (has_improved == false) {
//             	best_neighbour.normalise_ids();
//                 maxima.insert(best_neighbour);
//                 break;
//             }
//         }
//     }
// }

// }
