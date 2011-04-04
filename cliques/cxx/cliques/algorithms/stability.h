#ifndef CLIQUES_STABILITY_H
#define CLIQUES_STABILITY_H

#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <cliques/graphhelpers.h>

#include <cliques/structures/disjointset.h>

namespace cliques
{
/**
@brief  Functor for finding stability
*/
struct find_linearised_stability {
	std::vector<float> &markov_times_;

	find_linearised_stability(std::vector<float> &markov_times)
		:markov_times_(markov_times)
	{}

	template <typename T1, typename T2>
	void operator () (T1 &my_graph, T2 &my_partition,
				 std::vector<float> &stabilities)
	{
		float two_m = 2.0 * lemon::countEdges(my_graph);
		float first_term = 0.0;
		float second_term = 0.0;

		for (typename T2::PartIterator pitr = my_partition.begin(); pitr != my_partition.end(); ++pitr) {
			for (typename T2::NodeIterator n1itr = pitr.begin(); n1itr != pitr.end(); ++n1itr) {
				for (typename T2::NodeIterator n2itr = pitr.begin(); n2itr != pitr.end(); ++n2itr) {
					//std::cout << "n1: " << *n1itr << " n2: " << *n2itr << std::endl;
					int k1 = lemon::countIncEdges(my_graph, my_graph.nodeFromId(*n1itr));
					int k2 = lemon::countIncEdges(my_graph, my_graph.nodeFromId(*n2itr));
					float A = cliques::A(my_graph, *n1itr, *n2itr);
					first_term = first_term + float(k1*k2);
					second_term = second_term + A;
				}
			}
		}
		first_term = first_term/(two_m*two_m);
		second_term = second_term/two_m;

		float R;
		for (std::vector<float>::iterator t = markov_times_.begin(); t != markov_times_.end(); ++t)
		{
			R = 0.0;
			R = (1.0-*t) - first_term + *t * second_term;
			stabilities.push_back(R);
		}
	}

	template <typename G>
	void operator() (G &graph, cliques::DisjointSetForest<int> &partition,
			std::vector<float> &stabilities) {
		float two_m = 2.0 * lemon::countEdges(graph);
		float first_term = 0.0;
		float second_term = 0.0;

		for (typename G::NodeIt n1(graph); n1 != lemon::INVALID; ++n1) {
			for (typename G::NodeIt n2 = n1; ++n2 != lemon::INVALID;) {
				if (partition.find_set(graph.id(n1)) == partition.find_set(graph.id(n2))) {
					int k1 = lemon::countIncEdges(graph, n1);
					int k2 = lemon::countIncEdges(graph, n2);
					float A = cliques::A(graph, graph.id(n1), graph.id(n2));
					first_term = first_term + float(k1*k2);
					second_term = second_term + A;
				}
			}
		}

		/*for (typename cliques::DisjointSetForest<int>::PartIterator pitr = partition.begin(); pitr != partition.end(); ++pitr) {
			for (typename cliques::DisjointSetForest<int>::NodeIterator n1itr = pitr.begin(); n1itr != pitr.end(); ++n1itr) {
				for (typename cliques::DisjointSetForest<int>::NodeIterator n2itr = pitr.begin(); n2itr != pitr.end(); ++n2itr) {
					//std::cout << "n1: " << *n1itr << " n2: " << *n2itr << std::endl;
					int k1 = lemon::countIncEdges(graph, graph.nodeFromId(*n1itr));
					int k2 = lemon::countIncEdges(graph, graph.nodeFromId(*n2itr));
					float A = cliques::A(graph, *n1itr, *n2itr);
					first_term = first_term + float(k1*k2);
					second_term = second_term + A;
				}
			}
		}*/
		first_term = first_term/(two_m*two_m);
		second_term = second_term/two_m;

		float R;
		for (std::vector<float>::iterator t = markov_times_.begin(); t != markov_times_.end(); ++t)
		{
			R = 0.0;
			R = (1.0-*t) - first_term + *t * second_term;
			stabilities.push_back(R);
		}
	}
};



/**
@brief  Functor for finding stability of weighted graph
*/
template <typename W>
struct find_weighted_linearised_stability {
	std::vector<float> &markov_times_;
	W &weights;

	find_weighted_linearised_stability(std::vector<float> &markov_times,
			W &weights) : markov_times_(markov_times), weights(weights)
	{}


	template <typename G, typename T2>
	void operator () (G &graph, T2 &my_partition,
				 std::vector<float> &stabilities)
	{
		std::cout << "USING DA ORIGINAL";
		float two_m = 2.0 * cliques::find_total_weight(graph, weights);
		float first_term = 0.0;
		float second_term = 0.0;

		for (typename T2::PartIterator pitr = my_partition.begin(); pitr != my_partition.end(); ++pitr) {
			for (typename T2::NodeIterator n1itr = pitr.begin(); n1itr != pitr.end(); ++n1itr) {
				for (typename T2::NodeIterator n2itr = pitr.begin(); n2itr != pitr.end(); ++n2itr) {
					//std::cout << "n1: " << *n1itr << " n2: " << *n2itr << std::endl;
					float k1 = cliques::find_weighted_degree(graph, weights, graph.nodeFromId(*n1itr));
					float k2 = cliques::find_weighted_degree(graph, weights, graph.nodeFromId(*n2itr));
					typename G::Edge edge = lemon::findEdge(graph, graph.nodeFromId(*n1itr), graph.nodeFromId(*n2itr));
					float A = weights[edge];
					first_term = first_term + k1*k2;
					second_term = second_term + A;
				}
			}
		}
		first_term = first_term/(two_m*two_m);
		second_term = second_term/two_m;

		float R;
		for (std::vector<float>::iterator t = markov_times_.begin(); t != markov_times_.end(); ++t)
		{
			R = 0.0;
			R = (1.0-*t) - first_term + *t * second_term;
			stabilities.push_back(R);
		}
	}
};

struct find_combinatorial_stability {
	std::vector<float> &markov_times_;

	find_combinatorial_stability(std::vector<float> &markov_times)
		:markov_times_(markov_times)
	{}

	template <typename T1, typename T2>
	void operator () (T1 &my_graph, T2 &my_partition,
				 std::vector<float> &stabilities)
	{
		float two_m = 2.0 * lemon::countEdges(my_graph);
		float first_term = 0.0;
		float second_term = 0.0;

		for (typename T2::PartIterator pitr = my_partition.begin(); pitr != my_partition.end(); ++pitr) {
			for (typename T2::NodeIterator n1itr = pitr.begin(); n1itr != pitr.end(); ++n1itr) {
				for (typename T2::NodeIterator n2itr = pitr.begin(); n2itr != pitr.end(); ++n2itr) {
					//std::cout << "n1: " << *n1itr << " n2: " << *n2itr << std::endl;
					int k1 = lemon::countIncEdges(my_graph, my_graph.nodeFromId(*n1itr));
					int k2 = lemon::countIncEdges(my_graph, my_graph.nodeFromId(*n1itr));
					float A = cliques::A(my_graph, *n1itr, *n2itr);
					first_term = first_term + float(k1*k2);
					second_term = second_term + A;
				}
			}
		}
		first_term = first_term/(two_m*two_m);
		second_term = second_term/two_m;

		float R;
		for (std::vector<float>::iterator t = markov_times_.begin(); t != markov_times_.end(); ++t)
		{
			R = 0.0;
			R = (1.0-*t) - first_term + *t * second_term;
			stabilities.push_back(R);
		}
	}
};

/*void find_combinatorial_stability(std::map<int,int> &partition_map, std::vector<float> markov_times, std::vector<double> &stabilities)
{
	std::vector<std::vector<int > > p = map_to_vector(partition_map);
	double one_over_N = 1.0/(g->n);
	double one_over_N_squared = one_over_N * one_over_N;
	double mean_k = 0.0;

	for (int i =0; i<g->n; ++i)
	{
		mean_k = mean_k + g->degrees[i];
	}
	mean_k = mean_k / g->n;
	int time = 0;

	//calc once m = "sum of all weights"
	//degrees some of all rows or columns
	//Q = Q + A(p[s][i],p[s][j]) - float(g->degrees[p[s][i]]*g->degrees[p[s][j]])/(float(2*g->m));

	for (std::vector<float>::iterator t = markov_times.begin(); t != markov_times.end(); ++t)
	{
		double R = 0.0;
		for (unsigned int s = 0;  s < p.size(); ++s)
		{
			for (unsigned int i = 0; i< p[s].size(); ++i)
			{
				for (unsigned int j = 0; j < p[s].size(); ++j)
				{
					R += expmats[time][p[s][i]*g->n + p[s][j]] - one_over_N_squared;
				}
			}
		}
		//if (*t < 0.00801 ){
		//cout << *t << ": " << R << endl;
		//}
		stabilities.push_back(R);
		time++;
	}
}*/
}

#endif
