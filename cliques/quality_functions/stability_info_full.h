/* Copyright (c) Michael Schaub - michael.schaub09@imperial.ac.uk, 2010-2011 */
#pragma once
#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <lemon/maps.h>

#include <cliques/helpers/helpers.h>
#include <cliques/helpers/graphhelpers.h>
#include <cliques/quality_functions/internals/linearised_internals_info.h>
#include <cliques/helpers/math.h>



namespace clq {

// Create exponential mutual information graph for full information stability
template<typename G, typename W>
void create_mutual_information_graph_from_graph(G &graph, W &weights,
		double markov_time) {

  // some typedefs
  typedef typename G::EdgeIt EdgeIt;
  typedef typename G::Edge Edge;
  typedef typename G::NodeIt NodeIt;
  typedef typename G::Node Node;

  // get number of nodes and preallocate null_model_vec
  unsigned int num_nodes = lemon::countNodes(graph);
  std::vector<double> node_weighted_degree(num_nodes, 0);
  lemon::SmartGraph exp_graph;
  lemon::SmartGraph::EdgeMap<double> exp_graph_weights(exp_graph);
  std::vector<double> minus_t_D_inv_L(num_nodes * num_nodes, 0);
  double two_m = 0;

// get weights of each node (temp stored in null_model_vec) and total graph weight
 for (EdgeIt edge(graph); edge != lemon::INVALID; ++edge) {
      int node_u_id = graph.id(graph.u(edge));
      int node_v_id = graph.id(graph.v(edge));

      // weight of edge
      double weight = weights[edge];

      // if selfloop, only half of the weight has to be considered
      if (node_u_id == node_v_id) {
          weight = weight / 2;
      }

      // add weight to node weight and sum up total weight
      node_weighted_degree[node_u_id] += weight;
      node_weighted_degree[node_v_id] += weight;
      two_m += 2 * weight;

  }

  //initialise matrix and set diagonals to minus identity
  for (unsigned int i = 0; i < num_nodes * num_nodes; ++i) {
      if (i % (num_nodes + 1) == 0) {
          minus_t_D_inv_L[i] = -1;
      } else {
          minus_t_D_inv_L[i] = 0;
      }

  }
  // fill in the rest
  for (EdgeIt e(graph); e != lemon::INVALID; ++e) {
      Node u = graph.u(e);
      Node v = graph.v(e);
      int node_id_u = graph.id(u);
      int node_id_v = graph.id(v);
      double weight_uv = weights[e];
      // B_uv
      minus_t_D_inv_L[node_id_v + num_nodes * node_id_u] += weight_uv
              / node_weighted_degree[node_id_v];
      // B_vu
      minus_t_D_inv_L[node_id_u + num_nodes * node_id_v] += weight_uv
              / node_weighted_degree[node_id_u];
  }
  //clq::print_collection(minus_t_D_inv_L);
  // create graph structure
  //reserve memory space for number of nodes
  exp_graph.reserveNode(num_nodes);
  exp_graph.reserveEdge((num_nodes * (num_nodes + 1)) / 2);

  // add nodes
  for (unsigned int i = 0; i < num_nodes; ++i) {
      exp_graph.addNode();
  }
  // copy graph
  lemon::graphCopy(exp_graph, graph).run();
  // call expokit
  std::vector<double> exp_graph_vec = clq::exp(minus_t_D_inv_L,
          markov_time, num_nodes);



  // create new graph out of matrix
  for (unsigned int i = 0; i < num_nodes; ++i) {
      for (unsigned int j = i; j < num_nodes; ++j) {
          double weight = exp_graph_vec[num_nodes * i + j]
                  * (node_weighted_degree[j] / two_m);
          if (weight > 0) {
              weight = weight * log((exp_graph_vec[num_nodes * i + j]
                      * (node_weighted_degree[j] / two_m))
                      / ((node_weighted_degree[j] / two_m)
                              * (node_weighted_degree[i] / two_m))) / log(2);
              //clq::output(weight);
              lemon::SmartGraph::Edge edge = graph.addEdge(
                      graph.nodeFromId(i), graph.nodeFromId(j));
              weights[edge] = weight;
          }

      }
  }

}

} // end namespace
