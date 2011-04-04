#ifndef CLIQUES_MODULARITY_H
#define CLIQUES_MODULARITY_H

#include <vector>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>
#include <cliques/graphhelpers.h>

namespace cliques
{
/**
@brief  Functor for finding modularity
*/
struct find_modularity {
	template <typename G, typename T2>
	void operator () (G &graph, T2 &my_partition,
				 std::vector<float> &modularity)
	{
		float two_m = 2.0 * lemon::countEdges(graph);

		std::vector<float> body;
		body.push_back(0.0);

		for (typename T2::PartIterator pitr = my_partition.begin(); pitr != my_partition.end(); ++pitr) {
			for (typename T2::NodeIterator n1itr = pitr.begin(); n1itr != pitr.end(); ++n1itr) {
				for (typename T2::NodeIterator n2itr = pitr.begin(); n2itr != pitr.end(); ++n2itr) {
					int k1 = lemon::countIncEdges(graph, graph.nodeFromId(*n1itr));
					int k2 = lemon::countIncEdges(graph, graph.nodeFromId(*n1itr));
					float A = cliques::A(graph, *n1itr, *n2itr);
					float new_body =  A - float(k1*k2) / two_m;
                	//std::cout << "n1: " << *n1itr << " n2: " << *n2itr << "body: " << new_body << std::endl;

					body[0] = body[0] + A - float(k1*k2) / two_m;
				}
			}
		}
		body[0] = body[0] / two_m;
		modularity = body;
	}
};

}

#endif
