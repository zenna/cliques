/* Emulates partition_prune_fast.py */
/* Only works on a machine with (unsigned long) == 8 */
/* for graphs with at most 64 vertices */

// TODO
// Don't revisit nodes * done
// Remove obvious inefficencies
// Store basins
// Generalise to any graph

// Serialise Graph
// Extend to larger graph


// Stability
// Generalise code to compute for stability (and other measures) without refinding partitiosn/ refinding neighbours etc

// Random walk
// DG graph BSP
// Renaud method
//
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <bitset>
//#include <unordered_map>
//#include <unordered_set>
#include <set>
#include <map>
#include <boost/functional/hash.hpp>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>


#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#define MAX_LINE_LENGTH 1000
#define SET_MAX_SIZE 63

#define N 128
#define NUM_NODES 16
#define BITS_PER_NODE 4
using namespace std;

int num_sets = 0;
int num_partitions = 0;

//Global Vars
typedef std::vector<std::vector<int> > vec_2d_ints;
std::vector<unsigned long> partition_vals;
ofstream partitions_file;

void reorder_map(std::map<int,int> &partition_map) {
    int start_part = 0;
    std::map<int,int> temp_partition_map = partition_map;

    for (std::map<int,int>::iterator itr = temp_partition_map.begin(); itr != temp_partition_map.end(); ++itr) {
        int reference_part = itr->second;
	    std::vector<int> to_remove;

        for (std::map<int,int>::iterator new_itr = temp_partition_map.begin(); new_itr != temp_partition_map.end(); ++new_itr) {
            if (new_itr->second == reference_part) {
                partition_map[new_itr->first] = start_part;

                if (itr->first != new_itr->first)
                   to_remove.push_back(new_itr->first);
            }
        }
        for (int i =0; i< to_remove.size(); ++i) {
            temp_partition_map.erase(to_remove[i]);
        }
        start_part++;
    }
}

void print_2d_vector (vec_2d_ints best_partition) {
    vec_2d_ints::iterator itr;
    for ( itr = best_partition.begin(); itr != best_partition.end(); ++itr ) {
        std::vector<int>::iterator new_itr;
        for ( new_itr = itr->begin(); new_itr != itr->end(); ++new_itr ) {
            std::cout << (*new_itr)+1 << " ";
        }
        cout << "\n";
    }
}

void print_map(std::map<int,int> my_map, bool add_one = false) {

    for (std::map<int,int>::iterator itr = my_map.begin(); itr != my_map.end(); ++itr) {
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

struct Partition // the event log data structure
{
	int key;
	//std::vector<std::vector<char> > partitions;
	std::bitset<N> bit_partitions;
	float modularity;
};

struct bitset_hash
{
    size_t operator()(std::bitset<N> bitset) const
    {
        return boost::hash<std::string>()(bitset.to_string());
    }
};

typedef boost::unordered_map<std::bitset<N>, Partition*, bitset_hash> umap;

struct iterator_hash
{
    size_t operator()(umap::const_iterator it) const
    {
        return boost::hash<std::string>()(it->first.to_string());
    }
};

umap partition_map;

std::vector<Partition> partitions_vec;

struct SortByMod // comparison function
{
	bool operator() (
	const Partition & a,
	const Partition & b) const {
		return a.modularity < b.modularity;
	}

	static Partition min_value()
	{
	Partition dummy;
	dummy.modularity = -1.0;
	return dummy;
	}

	static Partition max_value()
	{
	Partition dummy;
	dummy.modularity = 1.0;
	return dummy;
	}
};

typedef struct graph
{
    int n;
    int m;
    int **links;
    int *degrees;
    int *capacities;
} graph;

/* size is the max number of sets */
/* current is the index of the first free set */
typedef struct zzpartition
{
    int size;
    int current;
    unsigned long *sets;
} zzpartition;

/* Global variables */
graph *g;
unsigned long nodes_to_place;
zzpartition *part;
unsigned long part_union;

void report_error(char *s)
{
    fprintf(stderr,"%s\n",s);
    exit(-1);
}

void free_graph(graph *g)
{
    if (g!=NULL)
    {
        if (g->links!=NULL)
        {
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

graph *graph_from_file(FILE *f)
{
    char line[MAX_LINE_LENGTH];
    int i, u, v;
    graph *g;

    assert( (g=(graph *)malloc(sizeof(graph))) != NULL );

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
    for(i=0; i<g->n; i++)
    {
        if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
            report_error("graph_from_file; read error (fgets) 2");
        if( sscanf(line, "%d %d\n", &v, &(g->capacities[i])) != 2 )
            report_error("graph_from_file; read error (sscanf) 2");
        if( v != i )
        {
            fprintf(stderr,"Line just read : %s\n i = %d; v = %d\n",line,i,v);
            report_error("graph_from_file: error while reading degrees");
        }
    }

    /* compute the number of links */
    g->m=0;
    for(i=0; i<g->n; i++)
        g->m += g->capacities[i];
    g->m /= 2;

    /* create contiguous space for links */
    if (g->n==0)
    {
        g->links = NULL;
        g->degrees = NULL;
        g->capacities = NULL;
    }
    else
    {
        if( (g->links=(int **)malloc(g->n*sizeof(int*))) == NULL )
            report_error("graph_from_file: malloc() error 3");
        if( (g->links[0]=(int *)malloc(2*g->m*sizeof(int))) == NULL )
            report_error("graph_from_file: malloc() error 4");
        for(i=1; i<g->n; i++)
            g->links[i] = g->links[i-1] + g->capacities[i-1];
    }

    /* read the links */
    for(i=0; i<g->m; i++)
    {
        if( fgets(line,MAX_LINE_LENGTH,f) == NULL )
            report_error("graph_from_file; read error (fgets) 3");
        if( sscanf(line, "%d %d\n", &u, &v) != 2 )
        {
            fprintf(stderr,"Attempt to scan link #%d failed. Line read:%s\n", i, line);
            report_error("graph_from_file; read error (sscanf) 3");
        }
        if ( (u>=g->n) || (v>=g->n) || (u<0) || (v<0) )
        {
            fprintf(stderr,"Line just read: %s",line);
            report_error("graph_from_file: bad node number");
        }
        if ( (g->degrees[u]>=g->capacities[u]) ||
                (g->degrees[v]>=g->capacities[v]) )
        {
            fprintf(stderr, "reading link %s\n", line);
            report_error("graph_from_file: too many links for a node");
        }
        g->links[u][g->degrees[u]] = v;
        g->degrees[u]++;
        g->links[v][g->degrees[v]] = u;
        g->degrees[v]++;
    }
    for(i=0; i<g->n; i++)
        if (g->degrees[i]!=g->capacities[i])
            report_error("graph_from_file: capacities <> degrees");
    /*  printf("%s\n", line);*/
    /* horrible hack */
    /*   if( fgets(line,MAX_LINE_LENGTH,f) != NULL ) */
    /*     printf("%s\n", line); */
    /*     report_error("graph_from_file; too many lines"); */

    return(g);
}

/* set functions */
unsigned long empty_set()
{
    return(0);
}

/* fails if e is already in set */
unsigned long add_to_set(int e, unsigned long s)
{
    unsigned long tmp=1;

    tmp=tmp<<e;
    assert( (s & tmp) == 0 );
    s += tmp;
    return(s);
}

unsigned long add_to_set_no_fail(int e, unsigned long s)
{
    unsigned long tmp=1;

    tmp=tmp<<e;
    if( (s & tmp) == 0 )
        return(s+tmp);
    return(s);
}

/* checks whether e is in the set */
/* fails if it is not */
unsigned long remove_from_set(int e, unsigned long s)
{
    unsigned long tmp=1;

    tmp=tmp<<e;
    assert( (s & tmp) != 0 );
    s -= tmp;
    return(s);
}

unsigned long set_union(unsigned long s1, unsigned long s2)
{
    return(s1 | s2);
}

/* returns a set with elements of s1 which are not in s2 */
unsigned long set_difference(unsigned long s1, unsigned long s2)
{
    return(s1 - (s1 & s2));
}

/* n is the largest possible element in the set */
void print_setold(unsigned long s, int n)
{
    int i;
    unsigned long test;

    test=1;
    i=0;
    while( i < n )
    {
        if( (test & s) != 0 )
            printf("%d ",i);
        i+=1;
        test=test<<1;
        /*      printf("i : %d, test : %ld\n", i, test);  */
    }
    if( (test & s) != 0 )
        printf("%d", i);
    printf("\n");
}

/* returns smallest element from non-empty set */
int set_get_first(unsigned long s)
{
    int i;
    unsigned long test;

    assert( s != 0 );
    test=1;
    i=0;
    while(1)
    {
        if( (test & s) != 0 )
            return(i);
        i+=1;
        test=test<<1;
    }
}

/* set has length one is corresponding number is a power of two */
/* http://www.cprogramming.com/tutorial/powtwosol.html */
int has_length_one(unsigned long x)
{
    return !((x-1) & x);
}

//Sets the value of a bitset to value at offset
void set_bitset(std::bitset<N> &bit_array, int &offset, int value) {
    std::bitset<BITS_PER_NODE> temp_bit_array = value;
    for (int i =0;i<BITS_PER_NODE;i++) {
        int j = offset + i;
        bit_array[j] = temp_bit_array[i];
    }
    offset = offset + BITS_PER_NODE;
}

void partition_long_to_map(unsigned long s, int n, int set_num, std::map<int, int> &temp_partition_map)
{
    int i;
    unsigned long test;

    test=1;
    i=0;
    while( i < n )
    {
        if( (test & s) != 0 ) {
            temp_partition_map[i] = set_num;
        }

        i+=1;
        test=test<<1;
        /*      printf("i : %d, test : %ld\n", i, test);  */
    }
    if( (test & s) != 0 ) {
        //printf("%d ", i);
        temp_partition_map[i] = set_num;
    }
}

zzpartition *create_partition(int size)
{
    zzpartition *p;

    p=(zzpartition *)malloc(sizeof(zzpartition));
    p->size=size;
    p->current=0;
    p->sets=(unsigned long *)calloc(size,sizeof(unsigned long));
    return(p);
}

zzpartition *add_set(unsigned long s, zzpartition *p)
{
    assert( p->current < p->size );
    p->sets[p->current]=s;
    p->current++;
    return(p);
}

zzpartition *remove_last_set(zzpartition *p)
{
    assert( p->current > 0 );
    p->sets[p->current]=empty_set();
    p->current--;
    return(p);
}

float A(int one, int two) {
    for(int i=0; i<g->degrees[one]; i++){
        if (g->links[one][i] == two) {
            return 1.0;
        }
    }

    return 0.0;
}

int get_bitset(std::bitset<N> &bit_array, int &offset) {
    std::bitset<BITS_PER_NODE> temp_bit_array;
    for (int i =0;i<BITS_PER_NODE;i++) {
        int j = offset + i;
        temp_bit_array[i] = bit_array[j];
    }
    offset = offset + BITS_PER_NODE;

    return int(temp_bit_array.to_ulong());

}

std::map<int,int> bitset_to_map(std::bitset<N> &bit_array) {
    int offset = 0;
    std::map<int,int> partition_map;

    for (int i = 0; i < NUM_NODES; i++) {
        int node = get_bitset(bit_array,offset);
        int part = get_bitset(bit_array,offset);
        partition_map[node] = part;
    }

    return partition_map;
};

std::bitset<N> map_to_bitset(std::map<int,int> &partition_map) {

    int offset = 0;
    std::bitset<N> bit_array;

    for (std::map<int, int>::iterator itr = partition_map.begin(); itr != partition_map.end(); ++itr) {
        set_bitset(bit_array, offset, itr->first);
        set_bitset(bit_array, offset, itr->second);
    }

    return bit_array;
}

std::vector<std::vector<int> > bitset_to_vector(std::bitset<N> &bit_array) {

    int offset = 0;
    std::vector<std::vector<int > > p;

    int first_part = 0;
    std::vector<int> temp_vec;
    for (int i = 0; i< NUM_NODES; i++) {
        int part_index = get_bitset(bit_array,offset);
        int node_index = get_bitset(bit_array,offset);
        if (part_index != first_part) {
            p.push_back(temp_vec);
            temp_vec.clear();
            first_part++;
        }
        temp_vec.push_back(node_index);
    }
    p.push_back(temp_vec);

    return p;
};

std::bitset<N> vector_to_bitset(vec_2d_ints &p) {

    std::bitset<N> temp_bitset;
    int offset = 0;
    for (int i = 0;i<p.size();++i) {
        for (int j = 0; j < p[i].size(); ++j) {
            set_bitset(temp_bitset,offset,i);
            set_bitset(temp_bitset,offset,p[i][j]);
        }
    }

    return temp_bitset;
};

std::set<std::set<int> > bitset_to_set(std::bitset<N> &bit_array) {

    int offset = 0;
    std::set<std::set<int > > p;

    int first_part = 0;
    std::set<int> temp_vec;
    for (int i = 0; i<16; i++) {
        int part_index = get_bitset(bit_array,offset);
        int node_index = get_bitset(bit_array,offset);
        if (part_index != first_part) {
            p.insert(temp_vec);
            temp_vec.clear();
            first_part++;
        }
        temp_vec.insert(node_index);
    }
    p.insert(temp_vec);

    return p;
}

vec_2d_ints map_to_vector(std::map<int,int> &partition_map) {

    int max_part_seen = -1;
    vec_2d_ints temp_partition_vector;
    for (std::map<int,int>::iterator itr = partition_map.begin(); itr != partition_map.end(); ++itr) {

        int part_num_diff = itr->second - max_part_seen;
        if (part_num_diff> 0) {
            for (int i=0;i<part_num_diff;++i) {
                std::vector<int> a;
                temp_partition_vector.push_back(a);
            }
            max_part_seen = itr->second;
        }
        temp_partition_vector[itr->second].push_back(itr->first);
    }

    return temp_partition_vector;
}

float find_modularity(std::map<int,int> &partition_map) {
    //Create std::vector from partition
    //int offset = 0;
    std::vector<std::vector<int > > p = map_to_vector(partition_map);

    /*int first_part = 0;
    std::vector<int> temp_vec;
    for (int i = 0; i<16; i++) {
        int part_index = get_bitset(bit_array,offset);
        int node_index = get_bitset(bit_array,offset);
        if (part_index != first_part) {
            p.push_back(temp_vec);
            temp_vec.clear();
            first_part++;
        }
        temp_vec.push_back(node_index);
    }
    p.push_back(temp_vec);*/

    float Q = 0.0;
    for (int s = 0;  s < p.size(); ++s) {
        for (int i = 0; i< p[s].size(); ++i) {
            for (int j = 0; j < p[s].size(); ++j) {
                Q = Q + A(p[s][i],p[s][j]) - float(g->degrees[p[s][i]]*g->degrees[p[s][j]])/(float(2*g->m));
            }
        }
    }
    Q = Q/(2*g->m);
    return Q;
}

void print_partition_old(zzpartition *p)
{
    int offset = 0;
    int i;
    std::map<int, int> temp_partition_map;
    Partition temp_partition;
    temp_partition.key = num_partitions;
    for(i=0; i<p->current; i++)
    {
        partition_long_to_map(p->sets[i],SET_MAX_SIZE,i, temp_partition_map);
        num_sets++;
        //printf("| ");
    }

    temp_partition.modularity = find_modularity(temp_partition_map);

    if (num_partitions % 100000 == 0) {
        cout << num_partitions << "\n";
    }

    //std::map<int, int> t_temp_partition_map = temp_partition_map ;

    reorder_map(temp_partition_map);

    for (std::map<int, int>::iterator itr = temp_partition_map.begin(); itr != temp_partition_map.end(); ++itr) {
        set_bitset(temp_partition.bit_partitions, offset, itr->first);
        set_bitset(temp_partition.bit_partitions, offset, itr->second);
    }

    /*std::string ref_string = "10001111001011100110110100101100011110110001101001101001001010000001011101010110010001010011010000100011000000100001000100000000";
    std::bitset<N> ref_bitset(ref_string);

    if (ref_bitset == temp_partition.bit_partitions) {
        std::cout << "waiting for youu!!\n";
        print_map(t_temp_partition_map);
        print_map(temp_partition_map);
        getchar();
    }*/

    partitions_vec.push_back(temp_partition);
    partition_map[partitions_vec.back().bit_partitions] = &partitions_vec.back();
    num_partitions++;
}

void print_partition(zzpartition *p)
{
    int i;
    num_partitions++;
    for(i=0; i<p->current; i++)
        num_sets++;
        // push_back vector of longs
        partitions_file << p->sets[i] << " ";
    partitions_file << "\n";
    if (num_partitions % 100000 == 0)
        cout << num_partitions << "\n";
}

/* generate connected partitions */

void gen_partitions_prune_aux(unsigned long comp_tmp,unsigned long neighbours,unsigned long forbidden)
{
    int v, vn,i;
    unsigned long new_neighbours, new_forbidden, new_comp_tmp;

    if( neighbours != 0 )
    {
        /* possibility to extend comp_tmp with a vertex from neighbours */

        vn = set_get_first(neighbours);
        comp_tmp=add_to_set(vn,comp_tmp);

        new_neighbours=remove_from_set(vn, neighbours);
        for(i=0; i<g->degrees[vn]; i++)
            new_neighbours=add_to_set_no_fail(g->links[vn][i],new_neighbours);
        new_neighbours=set_difference(new_neighbours,forbidden);

        forbidden=add_to_set(vn,forbidden);

        nodes_to_place=remove_from_set(vn,nodes_to_place);

        /*  Recursive call, extend comp_tmp with first neighbour" */
        gen_partitions_prune_aux(comp_tmp, new_neighbours, forbidden);

        nodes_to_place=add_to_set(vn,nodes_to_place);
        comp_tmp=remove_from_set(vn,comp_tmp);

        /* exploring the possibility to extend comp_tmp without using vn */
        new_neighbours=remove_from_set(vn,neighbours);

        /*  rec call, extend comp_tmp with other neighour */
        gen_partitions_prune_aux(comp_tmp, new_neighbours, forbidden);

        forbidden=remove_from_set(vn,forbidden);
    }
    else
    {
        /* all possibilities to extend comp_tmp have been considered, */
        /* new consider comp_tmp as a whole set */

        /* if comp_tmp has only one element we do not want it as a whole set */
        //if( has_length_one(comp_tmp) )
        //  return;

        if( nodes_to_place == 0 )
        {
            /* we have a zzpartition */
            part=add_set(comp_tmp,part);
            print_partition_old(part);
            remove_last_set(part);
            return;
        }
        /* start a new set with first vertex from nodes_to_place */

        v=set_get_first(nodes_to_place);
        nodes_to_place=remove_from_set(v,nodes_to_place);

        add_set(comp_tmp,part);
        part_union=set_union(comp_tmp,part_union);

        new_forbidden=add_to_set(v,part_union);
        new_comp_tmp=empty_set();
        new_comp_tmp=add_to_set(v,new_comp_tmp);

        new_neighbours=empty_set();
        for(i=0; i<g->degrees[v]; i++)
            new_neighbours=add_to_set_no_fail(g->links[v][i],new_neighbours);
        new_neighbours=set_difference(new_neighbours,new_forbidden);


        /* rec call, comp_tmp is a new set */
        gen_partitions_prune_aux(new_comp_tmp, new_neighbours, new_forbidden);

        nodes_to_place=add_to_set(v,nodes_to_place);

        part=remove_last_set(part);
        part_union=set_difference(part_union,comp_tmp);
    }
};

void print_1d_vector(std::vector<int> &vec) {

    for (std::vector<int>::iterator itr = vec.begin(); itr != vec.end(); ++itr ) {
        std::cout << *itr << " ";
    }

    std::cout << "\n";

}

bool is_connected(std::vector<int> part) {
    //Find neighbours only within partition
    //vec_2d_ints G;
    std::map<int, std::vector<int> > G;
    //std::cout << "part is =\n";
    //print_1d_vector(part);

    for (int i = 0; i<part.size(); ++i) {
        std::vector<int> temp_vector;
        for(int j=0; j< g->degrees[part[i]]; ++j) {
            if (std::find(part.begin(), part.end(), g->links[part[i]][j]) != part.end()) {
                temp_vector.push_back(g->links[part[i]][j]);
            }
         }
         G[part[i]] = temp_vector;
    }

    std::set<int> seen;
    std::set<int> nextlevel;
    nextlevel.insert(part[0]);  //# dict of nodes to check at next level
    while (nextlevel.size() != 0) {
        std::set<int> thislevel = nextlevel;  //# advance to next level
        nextlevel.clear();         //# and start a new list (fringe)
        for (std::set<int>::iterator v = thislevel.begin(); v != thislevel.end(); ++v){
            if (seen.find(*v) == seen.end()) {
                seen.insert(*v);// # set the level of vertex v
                for (int i = 0; i < G[*v].size(); ++i) {
                    nextlevel.insert(G[*v][i]); //# add neighbors of v
                }
            }
        }
    }
    return (seen.size() == part.size());  //# return all path lengths as dictionary
};

bool value_comparer(std::map<int,int>::value_type &i1, std::map<int,int>::value_type &i2) {
    return i1.second<i2.second;
}

umap::const_iterator find_neighbours(umap::const_iterator &centre, umap &partition_map, boost::unordered_set<umap::const_iterator, iterator_hash> &neighbours) {

    std::map<int,int> current_partition = bitset_to_map(centre->second->bit_partitions);
    float best_score = std::numeric_limits<int>::min();
    umap::const_iterator best_neighbour;

    // Highest of map values + 1 will be number of parts
    int num_parts = (std::max_element(current_partition.begin(), current_partition.end(),value_comparer))->second + 1;

    for (std::map<int,int>::iterator itr = current_partition.begin(); itr != current_partition.end(); ++itr) {
        for (int next_part = 0; next_part<num_parts+1; ++next_part) {
            std::map<int,int> neighbour_partition = current_partition;
            int dist_iterated = std::distance(current_partition.begin(), itr);
            std::map<int,int>::iterator n_itr = neighbour_partition.begin();
            std::advance (n_itr,dist_iterated);


            if (next_part != itr->second) {
                //std::cout << "next part is" << next_part << std::endl;
                int current_part = n_itr->second;
                n_itr->second = next_part;

                // Find only partitions changed in move set for following check of connectivity
                // More efficient than checking whole graph (I think?)
                std::vector<int> first_changed_part, last_changed_part;
                for (std::map<int,int>::iterator new_itr = neighbour_partition.begin(); new_itr != neighbour_partition.end(); ++new_itr) {
                    if (new_itr->second == current_part)
                        first_changed_part.push_back(new_itr->first);
                    if (new_itr->second == next_part)
                        last_changed_part.push_back(new_itr->first);
                }

                // Check the partition change created a connected graph
                if ( ((!first_changed_part.size()) || is_connected(first_changed_part) ) && ( (!last_changed_part.size()) || is_connected(last_changed_part) )) {
                    reorder_map(neighbour_partition);

                    std::bitset<N> temp_bit_array = map_to_bitset(current_partition);
                    std::bitset<N> new_bit_array = map_to_bitset(neighbour_partition);

                    if (partition_map[new_bit_array]->modularity > best_score) {
                        best_score = partition_map[new_bit_array]->modularity;
                        best_neighbour = partition_map.find(new_bit_array);
                        //std::cout << "best score so far is " << best_score << endl;
                    }

                    neighbours.insert(partition_map.find(new_bit_array));
                }
            }
        }
    }

    return best_neighbour;
}

using namespace lemon;
using namespace std;

typedef ListGraph::Node Node;
typedef ListGraph::Edge Edge;
typedef ListGraph::EdgeMap<int> LengthMap;
typedef lemon::ListGraph Graph;

//template <class float>
float find_bottleneck (ListGraph &g, CrossRefMap<ListGraph, Edge, float> &weight, Node &s, Node &t) {
    // Check if theres a single edge left
    if (countEdges(g) == 1) {
        CrossRefMap<ListGraph, Edge, float>::ValueIt itr = weight.beginValue();
        return *itr;
    }

    vector<int> deleted_edges;
    vector<Node> deleted_u;
    vector<Node> deleted_v;
    vector<float> deleted_weights;

    int num_edges = countEdges(g);
    int num_deleted = 0;

    // Assumes distinct weights?
    num_edges = countEdges(g);
    cout << "Edges before pruning: " << num_edges << endl;

    for (CrossRefMap<ListGraph, Edge, float>::ValueIt itr = weight.beginValue(); itr != weight.endValue(); ++itr) {
        if (num_deleted >= num_edges / 2) {
            break;
        }

        Edge temp_edge = weight(*itr);
        deleted_edges.push_back(g.id(temp_edge));
        deleted_u.push_back(g.u(temp_edge));
        deleted_v.push_back(g.v(temp_edge));
        deleted_weights.push_back(*itr);
        //g.erase(temp_edge);
        num_deleted++;
    }

    for (int i = 0;i<deleted_edges.size();++i) {
        cout << "deleteing" << deleted_edges[i] << endl;
        g.erase(g.edgeFromId(deleted_edges[i]));
    }

    num_edges = countEdges(g);
    cout << "Edges after pruning: " << num_edges << endl;

    getchar();

    Graph::NodeMap<int> comp(g);
    int num_components = connectedComponents(g,comp);

    if (comp[s] != comp[t]) {
        Graph temp_graph;
        CrossRefMap<ListGraph, Edge, float> temp_weights(temp_graph);
        for (int i=0;i<num_components;++i) {
            temp_graph.addNode();
        }
        s = temp_graph.nodeFromId(comp[s]);
        t = temp_graph.nodeFromId(comp[t]);

        for (unsigned int i = 0;i<deleted_v.size();++i) {
            int component_of_v = comp[deleted_v[i]];
            int component_of_u = comp[deleted_u[i]];
            if (component_of_u != component_of_v) {
                Node node_u = temp_graph.nodeFromId(component_of_u);
                Node node_v = temp_graph.nodeFromId(component_of_v);
                // Check if there is an edge
                Edge e = findEdge(temp_graph,node_u,node_v);
                if (e != INVALID) {
                    if (temp_weights[e] < deleted_weights[i]) {
                        temp_weights.set(e,deleted_weights[i]);
                    }
                }
                else {
                    Edge e = temp_graph.addEdge(node_u,node_v);
                    temp_weights.set(e,deleted_weights[i]);
                }
            }
        }
        return find_bottleneck(temp_graph,temp_weights,s,t);
    }
    else {
        return find_bottleneck(g,weight,s,t);
    }
}


int main(int argc, char* argv[]) {
    partitions_vec.reserve(3000000);
    unsigned long s,s2,comp_tmp,neighbours,forbidden;
    int i;
    zzpartition *p;

    //Initialisation
    FILE * pFile = fopen (argv[1],"r");
    g=graph_from_file(pFile);
    nodes_to_place=empty_set();
    for(i=1; i<g->n; i++)
        nodes_to_place=add_to_set(i,nodes_to_place);
    part=create_partition(g->n);
    part_union=empty_set();
    comp_tmp=empty_set();
    comp_tmp=add_to_set(0,comp_tmp);
    neighbours=empty_set();
    for(i=0; i<g->degrees[0]; i++)
        neighbours=add_to_set(g->links[0][i],neighbours);
    forbidden=empty_set();
    forbidden=add_to_set(0,forbidden);

    gen_partitions_prune_aux(comp_tmp, neighbours, forbidden);

    cout << "partitions: " << num_partitions << endl;

    boost::unordered_set<umap::const_iterator, iterator_hash> maxima;
    boost::unordered_set<umap::const_iterator, iterator_hash> visited;

    ListGraph g_landscape;
    for (i = 0;i<num_partitions;++i) {
        g_landscape.addNode();
    }
    CrossRefMap<ListGraph, Edge, float> cont_weights(g_landscape);

    int num_skipped = 0;
    int num_iterations = 0;
    for (umap::const_iterator part_itr = partition_map.begin(); part_itr != partition_map.end(); ++part_itr ) {
        // Check part_itr hasn't already been visited
        if (num_iterations % 1000 == 0) {
            std::cout << "the number of iterations is" << num_iterations << std::endl;
            std::cout << "so far number of maxima is " << maxima.size() << std::endl;
            std::cout << "size of visited is " << visited.size() << std::endl;
            std::cout << "num skipped " << num_skipped << std::endl;
        }

        // don't revisit places I've been to before (it's boring)
        if (visited.find(part_itr) != visited.end()) {
            num_skipped++;
            num_iterations++;
            continue;
        }
        umap::const_iterator p = part_itr;
        while (1) {
            boost::unordered_set<umap::const_iterator, iterator_hash> partition_neighbours;
            umap::const_iterator best_neighbour = find_neighbours(p,partition_map,partition_neighbours);
            cout << "node: " << p->second->key << " size: " << partition_neighbours.size() << endl; //From bitset find all neighbours, lookup in map, return iterator to best one
            visited.insert(p);

            // Create Graph of landscape from neighbours, edge_weight is lesser modularity of nodes
            Node node_u = g_landscape.nodeFromId(p->second->key);
            for (boost::unordered_set<umap::const_iterator, iterator_hash>::iterator itr = partition_neighbours.begin(); itr != partition_neighbours.end();++itr) {
                Node node_v = g_landscape.nodeFromId((*itr)->second->key);
                // Check if there is an edge
                Edge e = findEdge(g_landscape,node_u,node_v);
                if (e == INVALID) {
                    Edge e = g_landscape.addEdge(node_u,node_v);
                    float lesser_weight = p->second->modularity < (*itr)->second->modularity ? p->second->modularity : (*itr)->second->modularity;
                    cont_weights.set(e,lesser_weight);
                }
            }

            if (p->second->modularity >= best_neighbour->second->modularity) {
                maxima.insert(p);
                break;
            }
            else {
                if (visited.find(best_neighbour) != visited.end()) {
                    break;
                }
                p = best_neighbour;

            }
        }
        num_iterations++;
    }

    cout << "Number of Edges is " << countEdges(g_landscape) << " and is_conn: " << connected(g_landscape) << endl;
    std::cout << "Total number of maxima is " << maxima.size() << endl;

    for (boost::unordered_set<umap::const_iterator, iterator_hash>::iterator itr = maxima.begin(); itr !=  maxima.end(); ++itr) {
        bitset<N> temp_bitset = (*itr)->first;
        cout << (*itr)->second->modularity << " ";
        print_map(bitset_to_map(temp_bitset));
    }

    cout << "Initiating Graph Search" << countEdges(g_landscape) << " and is conn is " << connected(g_landscape) << endl;
    std::cout << "total number of maxima is " << maxima.size() << endl;

    for (boost::unordered_set<umap::const_iterator, iterator_hash>::iterator itr = maxima.begin(); itr !=  maxima.end(); ++itr) {
        for (boost::unordered_set<umap::const_iterator, iterator_hash>::iterator itr_2 = maxima.begin(); itr_2 !=  maxima.end(); ++itr_2) {
            if (itr != itr_2) {

                Graph g_landscape_temp;
                CrossRefMap<ListGraph, Edge, float> cont_weights_2(g_landscape_temp);
                graphCopy(g_landscape,g_landscape_temp).edgeMap(cont_weights,cont_weights_2).run();

                Node u = g_landscape.nodeFromId((*itr)->second->key);
                Node v = g_landscape.nodeFromId((*itr_2)->second->key);
                cout << find_bottleneck(g_landscape_temp,cont_weights_2,u,v) << "\t";
            }
            else {
                cout << 0.0 << "\t";
            }
        }
        cout <<endl;
    }

    free_graph(g);
    return 0;
}
