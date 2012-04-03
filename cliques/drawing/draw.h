#pragma once

#include <string>
#include <iostream>
#include <sstream>

#include <vector>
#include <map>
#include <gvc.h>

#include <cliques/structures/common.h>
#include <cliques/helpers.h>
#include <cliques/algorithms/scalers.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


namespace cliques
{

// TODO - set good defaults and provide mechanism to change these if desired
// provide way to input position data

/**
@brief  Draws the graph using graphviz
*/
class draw_graph{
private:
	Agraph_t *g;
	GVC_t *gvc;
	std::map<int,Agedge_t *> edge_to_gv_edge;
	cliques::linear_scaler scaler;

public:
	template <typename G>
	draw_graph(G &graph) {
		gvc = gvContext();
		g = agopen(&(cliques::string_to_char_p("g"))[0], AGRAPHSTRICT);
		for (typename G::NodeIt n(graph); n != lemon::INVALID; ++n) {
			std::string node_name = cliques::to_string(graph.id(n));
			agnode(g,&(cliques::string_to_char_p(node_name))[0]);
		}
		srand ( time(NULL) );

		for (typename G::EdgeIt itr(graph); itr != lemon::INVALID; ++itr) {
	        if (graph.id(itr) % 1000 == 0) {
	            std::cout << "edge" << graph.id(itr) << "\n";
	        }
			std::string n_name = cliques::to_string(graph.id(graph.u(itr)));
			std::string m_name = cliques::to_string(graph.id(graph.v(itr)));

			Agnode_t *n = agnode(g, &(cliques::string_to_char_p(n_name))[0]);
			Agnode_t *m = agnode(g, &(cliques::string_to_char_p(m_name))[0]);
			edge_to_gv_edge[graph.id(itr)] = agedge(g, n, m);

			/*float r = (float)rand()/(float)RAND_MAX;
			std::string rand = cliques::to_string(r);
			char* rand_b = new char[rand.length() + 1];
			strcpy(rand_b, rand.c_str());
			agsafeset(edge_to_gv_edge[graph.id(itr)], "len", rand_b, "");

			delete [] rand_b;*/

			scaler = cliques::linear_scaler(0.0,0.0);
		}
		set_defaults();
	}

	~draw_graph() {
		/* Free layout data */
		gvFreeLayout(gvc, g);
		/* Free graph structures */
		agclose(g);
		/* close output file, free context, and return number of errors */
		gvFreeContext(gvc);
	}

	void draw(std::string filename, std::string algorithm = "sfdp") {
		gvLayout (gvc, g, &(cliques::string_to_char_p(algorithm))[0]);
		gvRenderFilename (gvc, g, "png", &(cliques::string_to_char_p(filename))[0]);
	}

	template <typename CL>
	void add_node_map(CL make_node_colour_map) {
	    add_node_colours(make_node_colour_map);
	}

	template <typename CL>
	void add_edge_map(CL make_edge_colour_map) {
	    add_edge_colours(make_edge_colour_map);
	}

    /**
    @brief  Adds graph weights to internal graph so they can be used to draw

    @param[in]  graph     Reference to graph to be drawn
    @param[in]  weights Edgemap relating edge -> weight
    */
	template <typename G, typename W>
	void add_weights(G &graph, W &weights) {

		// Scale edge weights to 0-1 scale
		for (std::map<int,Agedge_t *>::iterator itr = edge_to_gv_edge.begin();
						itr != edge_to_gv_edge.end(); ++itr) {
			typename G::Edge edge = graph.edgeFromId(itr->first);
			scaler.update_maxima(float(weights[edge]));
		}

		for (std::map<int,Agedge_t *>::iterator itr = edge_to_gv_edge.begin();
				itr != edge_to_gv_edge.end(); ++itr) {

			typename G::Edge edge = graph.edgeFromId(itr->first);

			std::string weight_str = cliques::to_string(1/scaler(weights[edge]));
			agsafeset(itr->second, &(cliques::string_to_char_p("len"))[0], &(cliques::string_to_char_p(weight_str))[0], &(cliques::string_to_char_p(""))[0]);
		}
	}

private:
	template <typename CL>
	void add_node_colours (CL &make_node_colour_map) {
		// Add node colour attributes
		std::map<int, cliques::xyz_colour<float> > colour_map = make_node_colour_map();
	    for (std::map<int, cliques::xyz_colour<float> >::iterator itr = colour_map.begin();
	    		itr != colour_map.end(); ++itr) {

	    	std::string n_name = cliques::to_string(itr->first);
	    	Agnode_t *n = agnode(g, &(cliques::string_to_char_p(n_name))[0]);

	    	std::ostringstream hsv;  // Declare an output string stream.
	    	hsv << itr->second.x << " "
	    		<< itr->second.y << " "
	    		<< itr->second.z;

	    	agsafeset(n, &(cliques::string_to_char_p("color"))[0], &(cliques::string_to_char_p(hsv.str()))[0], &(cliques::string_to_char_p(""))[0]);
	    }
	}

	template <typename CL>
	void add_edge_colours (CL &make_edge_colour_map) {
		std::map<int, cliques::xyz_colour<float> > colour_map = make_edge_colour_map();
	    for (std::map<int, cliques::xyz_colour<float> >::iterator itr = colour_map.begin();
	    		itr != colour_map.end(); ++itr) {

	    	std::ostringstream hsv;  // Declare an output string stream.
	    	hsv << itr->second.x << " "
	    		<< itr->second.y << " "
	    		<< itr->second.z;

	    	agsafeset(edge_to_gv_edge[itr->first], &(cliques::string_to_char_p("color"))[0], &(cliques::string_to_char_p(hsv.str()))[0], &(cliques::string_to_char_p(""))[0]);
	    }
	}

    /**
    @brief  Set default graph drawing attributes

    Makes decent looking graphs generally without config
    */
	void set_defaults() {
		agraphattr(g, &(cliques::string_to_char_p("dpi"))[0], &(cliques::string_to_char_p("300"))[0]);
		agnodeattr(g, &(cliques::string_to_char_p("width"))[0], &(cliques::string_to_char_p(".05"))[0]);
		agnodeattr(g, &(cliques::string_to_char_p("height"))[0], &(cliques::string_to_char_p(".05"))[0]);
		agnodeattr(g, &(cliques::string_to_char_p("label"))[0], &(cliques::string_to_char_p(""))[0]);
		agedgeattr(g, &(cliques::string_to_char_p("penwidth"))[0], &(cliques::string_to_char_p("0.1"))[0]);
	}

};

}
