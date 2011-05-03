/* Copyright (c) Modified by Zenna Tavares-zennatavares@gmail.com, 2010-2011 */
#ifndef CLIQUES_ALLPARTITIONS_H
#define CLIQUES_ALLPARTITIONS_H
#include <assert.h>
#include <cstdlib>
#include <iostream>

#include <map>

#include <lemon/list_graph.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#include <cliques/structures/common.h>
#include <cliques/structures/partition.h>
#include <cliques/helpers.h>

namespace cliques {

typedef struct c_partition {
	int size;
	int current;
	unsigned long *sets;
} c_partition;

typedef struct c_graph {
	int n;
	int m;
	int **links;
	int *degrees;
	int *capacities;
} c_graph;

int set_get_first(unsigned long s) {
	int i;
	unsigned long test;

	assert( s != 0 );
	test = 1;
	i = 0;
	while (1) {
		if ((test & s) != 0)
			return (i);
		i += 1;
		test = test << 1;
	}
}

unsigned long add_to_set_no_fail(int e, unsigned long s) {
	unsigned long tmp = 1;

	tmp = tmp << e;
	if ((s & tmp) == 0)
		return (s + tmp);
	return (s);
}

/* checks whether e is in the set */
/* fails if it is not */
unsigned long remove_from_set(int e, unsigned long s) {
	unsigned long tmp = 1;

	tmp = tmp << e;
	assert( (s & tmp) != 0 );
	s -= tmp;
	return (s);
}

unsigned long set_union(unsigned long s1, unsigned long s2) {
	return (s1 | s2);
}

/* returns a set with elements of s1 which are not in s2 */
unsigned long set_difference(unsigned long s1, unsigned long s2) {
	return (s1 - (s1 & s2));
}

unsigned long empty_set() {
	return (0);
}

unsigned long add_to_set(int e, unsigned long s) {
	unsigned long tmp = 1;

	tmp = tmp << e;
	assert( (s & tmp) == 0 );
	s += tmp;
	return (s);
}

c_partition *create_partition(int size) {
	c_partition *p;

	p = (c_partition *) malloc(sizeof(c_partition));
	p->size = size;
	p->current = 0;
	p->sets = (unsigned long *) calloc(size + 1, sizeof(unsigned long));
	return (p);
}

c_partition *add_set(unsigned long s, c_partition *p) {
	assert( p->current < p->size );
	p->sets[p->current] = s;
	p->current++;
	return (p);
}

c_partition *remove_last_set(c_partition *p) {
	assert( p->current > 0 );
	//std::cout << "p->current =" << p->current << " size of sets is: " << p->size << endl;
	p->sets[p->current] = empty_set();
	p->current--;
	return (p);
}

void partition_long_to_map(unsigned long s, int n, int set_num,
		std::map<int, int> &temp_partition_map) {
	int i;
	unsigned long test;

	test = 1;
	i = 0;
	while (i < n) {
		if ((test & s) != 0) {
			temp_partition_map[i] = set_num;
		}

		i += 1;
		test = test << 1;
		/*      printf("i : %d, test : %ld\n", i, test);  */
	}
	if ((test & s) != 0) {
		//printf("%d ", i);
		temp_partition_map[i] = set_num;
	}
}

template<typename P>
void add_partition_to_set(
		c_partition *p,
		int num_nodes,
		int &num_partitions,
		boost::unordered_set<P, cliques::partition_hash,
				cliques::partition_equal> &all_partitions) {

	P new_partition(num_nodes);

	for (int j = 0; j < p->current; j++) {
		int s = p->sets[j];
		int n = SET_MAX_SIZE;
		int set_num = j;
		int i;
		unsigned long test;

		test = 1;
		i = 0;
		while (i < n) {
			if ((test & s) != 0) {
				new_partition.add_node_to_set(i, set_num);
				//temp_partition_map[i] = set_num;
			}

			i += 1;
			test = test << 1;
			/*      printf("i : %d, test : %ld\n", i, test);  */
		}
		if ((test & s) != 0) {
			//printf("%d ", i);
			new_partition.add_node_to_set(i, set_num);
			//temp_partition_map[i] = set_num;
		}
	}
	new_partition.normalise_ids();
	all_partitions.insert(new_partition);

	num_partitions++;
    if (num_partitions % 100 == 0) {
        std::cout << num_partitions << "\n";
    }
}

void save_partition(c_partition *p, cliques::umap &partition_map,
		int &num_partitions) {
	int offset = 0;
	int i;
	std::map<int, int> temp_partition_map;
	cliques::Partition temp_partition;
	temp_partition.key = num_partitions;
	for (i = 0; i < p->current; i++) {
		partition_long_to_map(p->sets[i], SET_MAX_SIZE, i, temp_partition_map);
	}
	//temp_partition.modularity = find_modularity(temp_partition_map);
	//find_stability_combinatorial(temp_partition_map,temp_partition.stability,m_times);
	//std::map<int, int> t_temp_partition_map = temp_partition_map;
	cliques::reorder_map(temp_partition_map);
	for (std::map<int, int>::iterator itr = temp_partition_map.begin(); itr
			!= temp_partition_map.end(); ++itr) {
		cliques::set_bitset(temp_partition.bit_partitions, offset, itr->first);
		cliques::set_bitset(temp_partition.bit_partitions, offset, itr->second);
	}

	//partitions_vec.push_back(temp_partition);
	partition_map[temp_partition.bit_partitions] = temp_partition;// &partitions_vec.back();
	num_partitions++;
	if (num_partitions % 100000 == 0) {
		std::cout << num_partitions << "\n";
	}
}

template<typename P>
void gen_partitions_prune_aux(
		unsigned long comp_tmp,
		unsigned long neighbours,
		unsigned long forbidden,
		c_graph* &g,
		c_partition*& part,
		unsigned long & nodes_to_place,
		unsigned long & part_union,
		cliques::umap & partition_map,
		int & num_partitions,
		boost::unordered_set<P, cliques::partition_hash,
				cliques::partition_equal> &all_partitions) {
	int v, vn, i;
	unsigned long new_neighbours, new_forbidden, new_comp_tmp;

	if (neighbours != 0) {
		/* possibility to extend comp_tmp with a vertex from neighbours */

		vn = set_get_first(neighbours);
		comp_tmp = add_to_set(vn, comp_tmp);

		new_neighbours = remove_from_set(vn, neighbours);
		for (i = 0; i < g->degrees[vn]; i++)
			new_neighbours
					= add_to_set_no_fail(g->links[vn][i], new_neighbours);
		new_neighbours = set_difference(new_neighbours, forbidden);

		forbidden = add_to_set(vn, forbidden);

		nodes_to_place = remove_from_set(vn, nodes_to_place);

		/*  Recursive call, extend comp_tmp with first neighbour" */
		gen_partitions_prune_aux(comp_tmp, new_neighbours, forbidden, g, part,
				nodes_to_place, part_union, partition_map, num_partitions,
				all_partitions);

		nodes_to_place = add_to_set(vn, nodes_to_place);
		comp_tmp = remove_from_set(vn, comp_tmp);

		/* exploring the possibility to extend comp_tmp without using vn */
		new_neighbours = remove_from_set(vn, neighbours);

		/*  rec call, extend comp_tmp with other neighour */
		gen_partitions_prune_aux(comp_tmp, new_neighbours, forbidden, g, part,
				nodes_to_place, part_union, partition_map, num_partitions,
				all_partitions);
		forbidden = remove_from_set(vn, forbidden);
	} else {
		/* all possibilities to extend comp_tmp have been considered, */
		/* new consider comp_tmp as a whole set */

		/* if comp_tmp has only one element we do not want it as a whole set */
		//if( has_length_one(comp_tmp) )
		//  return;

		if (nodes_to_place == 0) {
			/* we have a c_partition */
			part = add_set(comp_tmp, part);
			//cout << "ERROR " << num_partitions << endl;
			//print_partition_old(part);
			//save_partition(part, partition_map, num_partitions);
			add_partition_to_set(part, g->n, num_partitions, all_partitions);
			//std::cout << "found\n";
			remove_last_set(part);
			return;
		}
		/* start a new set with first vertex from nodes_to_place */

		v = set_get_first(nodes_to_place);
		nodes_to_place = remove_from_set(v, nodes_to_place);

		add_set(comp_tmp, part);
		part_union = set_union(comp_tmp, part_union);

		new_forbidden = add_to_set(v, part_union);
		new_comp_tmp = empty_set();
		new_comp_tmp = add_to_set(v, new_comp_tmp);

		new_neighbours = empty_set();
		for (i = 0; i < g->degrees[v]; i++)
			new_neighbours = add_to_set_no_fail(g->links[v][i], new_neighbours);
		new_neighbours = set_difference(new_neighbours, new_forbidden);

		/* rec call, comp_tmp is a new set */
		gen_partitions_prune_aux(new_comp_tmp, new_neighbours, new_forbidden,
				g, part, nodes_to_place, part_union, partition_map,
				num_partitions, all_partitions);

		nodes_to_place = add_to_set(v, nodes_to_place);

		part = remove_last_set(part);
		part_union = set_difference(part_union, comp_tmp);
	}
}

void report_error(char *s){
  fprintf(stderr,"%s\n",s);
  exit(-1);
}

c_graph *graph_from_file(FILE *f){
  char line[MAX_LINE_LENGTH];
  int i, u, v;
  c_graph *g;

  assert( (g=(c_graph *)malloc(sizeof(c_graph))) != NULL );

  /* read n */
  if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
    report_error("graph_from_file: read error (fgets) 1");
  if( sscanf(line, "%d\n", &(g->n)) != 1 )
    report_error("graph_from_file: read error (sscanf) 2");

  /* read the degree sequence */
  if( (g->capacities=(int *)malloc(g->n*sizeof(int))) == NULL )
    report_error("graph_from_file: malloc() error 2");
  if( (g->degrees=(int *)calloc(g->n,sizeof(int))) == NULL )
    report_error("graph_from_file: calloc() error");
  for(i=0;i<g->n;i++){
    if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
      report_error("graph_from_file; read error (fgets) 2");
    if( sscanf(line, "%d %d\n", &v, &(g->capacities[i])) != 2 )
      report_error("graph_from_file; read error (sscanf) 2");
    if( v != i ){
      fprintf(stderr,"Line just read : %s\n i = %d; v = %d\n",line,i,v);
      report_error("graph_from_file: error while reading degrees");
    }
  }

  /* compute the number of links */
  g->m=0;
  for(i=0;i<g->n;i++)
    g->m += g->capacities[i];
  g->m /= 2;

  /* create contiguous space for links */
  if (g->n==0){
    g->links = NULL; g->degrees = NULL; g->capacities = NULL;
  }
  else {
    if( (g->links=(int **)malloc(g->n*sizeof(int*))) == NULL )
      report_error("graph_from_file: malloc() error 3");
    if( (g->links[0]=(int *)malloc(2*g->m*sizeof(int))) == NULL )
      report_error("graph_from_file: malloc() error 4");
    for(i=1;i<g->n;i++)
      g->links[i] = g->links[i-1] + g->capacities[i-1];
  }

  /* read the links */
  for(i=0;i<g->m;i++) {
    if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
      report_error("graph_from_file; read error (fgets) 3");
    if( sscanf(line, "%d %d\n", &u, &v) != 2 ){
      fprintf(stderr,"Attempt to scan link #%d failed. Line read:%s\n", i, line);
      report_error("graph_from_file; read error (sscanf) 3");
    }
    if ( (u>=g->n) || (v>=g->n) || (u<0) || (v<0) ) {
      fprintf(stderr,"Line just read: %s",line);
      report_error("graph_from_file: bad node number");
    }
    if ( (g->degrees[u]>=g->capacities[u]) ||
         (g->degrees[v]>=g->capacities[v]) ){
      fprintf(stderr, "reading link %s\n", line);
      report_error("graph_from_file: too many links for a node");
    }
    g->links[u][g->degrees[u]] = v;
    g->degrees[u]++;
    g->links[v][g->degrees[v]] = u;
    g->degrees[v]++;
  }
  for(i=0;i<g->n;i++)
    if (g->degrees[i]!=g->capacities[i])
      report_error("graph_from_file: capacities <> degrees");
  /*  printf("%s\n", line);*/
  /* horrible hack */
/*   if( fgets(line,MAX_LINE_LENGTH,f) != NULL ) */
/*     printf("%s\n", line); */
/*     report_error("graph_from_file; too many lines"); */

  return(g);
}

void free_graph(c_graph *g){
  if (g!=NULL) {
    if (g->links!=NULL) {
      if (g->links[0]!=NULL)
        free(g->links[0]);
      free(g->links);
    }
    if (g->capacities!=NULL)
      free(g->capacities);
    if (g->degrees!=NULL)
      free(g->degrees);
    free(g);
  }
}

/**
 @brief  Find all connected partitions of a graph

 TODO remove mem leak (delete) from c_partition
 TODO save partitions into partition map
 */
template<typename G, typename P>
void find_connected_partitions(
		G &graph,
		boost::unordered_set<P, cliques::partition_hash,
				cliques::partition_equal> &all_partitions) {
	cliques::umap partition_map;
	/*c_graph *g = new c_graph;
	g->n = lemon::countNodes(graph);
	g->m = lemon::countEdges(graph);
	g->degrees = new int[g->n];
	g->links = new int *[g->n]; // alloc(g->n*sizeof(int*))) == NULL )

	// Convert Lemon Graph to c_graph
	for (typename G::NodeIt itr(graph); itr != lemon::INVALID; ++itr) {
		int degree = 0;
		for (typename G::IncEdgeIt e_itr(graph, itr); e_itr != lemon::INVALID; ++e_itr) {
			degree++;
		}
		g->degrees[graph.id(itr)] = degree;

		g->links[graph.id(itr)] = new int[degree];

		int i = 0;
		for (typename G::IncEdgeIt e_itr(graph, itr); e_itr != lemon::INVALID; ++e_itr) {
			g->links[graph.id(itr)][i] = graph.id(graph.runningNode(e_itr));
			++i;
		}
	}*/

	c_graph *g;

	unsigned long nodes_to_place;
	c_partition *part;
	unsigned long part_union;
	unsigned long comp_tmp, neighbours, forbidden;
	int num_partitions = 0;

	FILE *pFile;
	pFile=fopen("/home/zenna/repos/graph-codes/cliques/data/graphs/barbell_n12.nke", "r");
	g=graph_from_file(pFile);

	nodes_to_place = empty_set();
	for (int i = 1; i < g->n; i++)
		nodes_to_place = add_to_set(i, nodes_to_place);
	part = create_partition(g->n);
	part_union = empty_set();

	comp_tmp = empty_set();
	comp_tmp = add_to_set(0, comp_tmp);
	neighbours = empty_set();
	for (int i = 0; i < g->degrees[0]; i++)
		neighbours = add_to_set(g->links[0][i], neighbours);
	forbidden = empty_set();
	forbidden = add_to_set(0, forbidden);

	gen_partitions_prune_aux(comp_tmp, neighbours, forbidden, g, part,
			nodes_to_place, part_union, partition_map, num_partitions,
			all_partitions);

	free_graph(g);
	/*delete[] g->degrees;
	for (int i = 0; i < g->n; ++i) {
		delete[] g->links[i];
	}
	delete[] g->links;
	delete g;*/
}

}

#endif
