/* Emulates partition_prune_fast.py */
/* Only works on a machine with (unsigned long) == 8 */
/* for graphs with at most 64 vertices */
#include <cstdlib>
#include <lemon/lgf_writer.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <bitset>
#include <set>
#include <map>
#include <math.h>


#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

#if defined USE_BOOST
    #include <boost/functional/hash.hpp>
    #include <boost/unordered_map.hpp>
    #include <boost/unordered_set.hpp>
#else
    #include <unordered_map>
    #include <unordered_set>
#endif

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
        for (unsigned int i =0; i< to_remove.size(); ++i) {
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
	double modularity;
	std::vector<double> stability;
};

#if defined USE_BOOST
typedef boost::hash<std::string> shash;
#else
typedef std::hash<std::string> shash;
#endif

struct bitset_hash
{
    size_t operator()(std::bitset<N> bitset) const
    {
        return shash()(bitset.to_string());
    }
};

#if defined USE_BOOST
typedef boost::unordered_map<std::bitset<N>, Partition*, bitset_hash> umap;
#else
typedef std::unordered_map<std::bitset<N>, Partition*, bitset_hash> umap;
#endif

struct iterator_hash
{
    size_t operator()(umap::const_iterator it) const
    {
        return shash()(it->first.to_string());
    }
};

#if defined USE_BOOST
typedef boost::unordered_set<umap::const_iterator, iterator_hash> uset;
#else
typedef std::unordered_set<umap::const_iterator, iterator_hash> uset;
#endif

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
    p->sets=(unsigned long *)calloc(size+1,sizeof(unsigned long));
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
    //std::cout << "p->current =" << p->current << " size of sets is: " << p->size << endl;
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

float K(int one, int two) {
    // if connected return degree of j

    for(int i=0; i<g->degrees[one]; i++){
        if (g->links[one][i] == two) {
            return g->degrees[two];
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
    for (unsigned int i = 0;i<p.size();++i) {
        for (unsigned int j = 0; j < p[i].size(); ++j) {
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
    for (unsigned int s = 0;  s < p.size(); ++s) {
        for (unsigned int i = 0; i< p[s].size(); ++i) {
            for (unsigned int j = 0; j < p[s].size(); ++j) {
                Q = Q + A(p[s][i],p[s][j]) - float(g->degrees[p[s][i]]*g->degrees[p[s][j]])/(float(2*g->m));
            }
        }
    }
    Q = Q/(2*g->m);
    return Q;
}

void find_stability(std::map<int,int> &partition_map, std::vector<float> &stabilities, std::vector<float> markov_times) {
    std::vector<std::vector<int > > p = map_to_vector(partition_map);
    float two_m = 2*g->m;
    float first_term = 0.0;
    float second_term = 0.0;

    for (unsigned int s = 0;  s < p.size(); ++s) {
        for (unsigned int i = 0; i< p[s].size(); ++i) {
            for (unsigned int j = 0; j < p[s].size(); ++j) {
                first_term = first_term + float(g->degrees[p[s][i]]*g->degrees[p[s][j]]);
                second_term = second_term + A(p[s][i],p[s][j]);
            }
        }
    }

    first_term = first_term/(two_m*two_m);
    second_term = second_term/two_m;

    float R;
    for (std::vector<float>::iterator t = markov_times.begin(); t != markov_times.end(); ++t) {
        R = 0.0;
        R = (1-*t) - first_term + *t * second_term;
        stabilities.push_back(R);
    }
};

std::vector<std::vector<double> > expmats;
std::vector<double> m_comb;
std::vector<std::vector<double> > k_comb;

void find_stability_combinatorial(std::map<int,int> &partition_map, std::vector<double> &stabilities, std::vector<float> markov_times) {
    std::vector<std::vector<int > > p = map_to_vector(partition_map);
    double one_over_N = 1.0/(g->n);
    double one_over_N_squared = one_over_N * one_over_N;
    double mean_k = 0.0;

    for (int i =0;i<g->n;++i) {
        mean_k = mean_k + g->degrees[i];
    }
    mean_k = mean_k / g->n;



    int time = 0;

    //calc once m = "sum of all weights"
    //degrees some of all rows or columns


    //Q = Q + A(p[s][i],p[s][j]) - float(g->degrees[p[s][i]]*g->degrees[p[s][j]])/(float(2*g->m));


    for (std::vector<float>::iterator t = markov_times.begin(); t != markov_times.end(); ++t) {
        double R = 0.0;
        for (unsigned int s = 0;  s < p.size(); ++s) {
            for (unsigned int i = 0; i< p[s].size(); ++i) {
                for (unsigned int j = 0; j < p[s].size(); ++j) {
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
}

std::vector<float> m_times;

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
    find_stability_combinatorial(temp_partition_map,temp_partition.stability,m_times);

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
            //cout << "ERROR " << num_partitions << endl;
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

    for (unsigned int i = 0; i<part.size(); ++i) {
        std::vector<int> temp_vector;
        for(int j=0; j< g->degrees[part[i]]; ++j) {
            if (std::find(part.begin(), part.end(), g->links[part[i]][j]) != part.end()) {
                temp_vector.push_back(g->links[part[i]][j]);
            }
         }
         G[part[i]] = temp_vector;
    }

    /*std::cout << "\nG is:\n";
    for (std::map<int, std::vector<int> >::iterator itr = G.begin(); itr != G.end(); ++itr ) {
        std::cout << itr->first << ": ";
        print_1d_vector(itr->second);
    }*/

    //print_2d_vector(G);
    std::set<int> seen;
    std::set<int> nextlevel;
    nextlevel.insert(part[0]);  //# dict of nodes to check at next level
    while (nextlevel.size() != 0) {
        std::set<int> thislevel = nextlevel;  //# advance to next level
        nextlevel.clear();         //# and start a new list (fringe)
        for (std::set<int>::iterator v = thislevel.begin(); v != thislevel.end(); ++v){
            if (seen.find(*v) == seen.end()) {
                seen.insert(*v);// # set the level of vertex v
                for (unsigned int i = 0; i < G[*v].size(); ++i) {
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

void find_neighbours(umap::const_iterator &centre, umap &partition_map, uset &neighbours) {

    std::map<int,int> current_partition = bitset_to_map(centre->second->bit_partitions);
    //float best_score = std::numeric_limits<int>::min();
    //umap::const_iterator best_neighbour;

    // Highest of map values + 1 will be number of parts
    int num_parts = (std::max_element(current_partition.begin(), current_partition.end(),value_comparer))->second + 1;

    for (std::map<int,int>::iterator itr = current_partition.begin(); itr != current_partition.end(); ++itr) {
        for (int next_part = 0; next_part<num_parts+1; ++next_part) {
            std::map<int,int> neighbour_partition = current_partition;
            int dist_iterated = std::distance(current_partition.begin(), itr);
            std::map<int,int>::iterator n_itr = neighbour_partition.begin();
            std::advance (n_itr,dist_iterated);

            if (next_part != itr->second) {
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
                    std::bitset<N> new_bit_array = map_to_bitset(neighbour_partition);
                    neighbours.insert(partition_map.find(new_bit_array));
                }
            }
        }
    }
}

using namespace lemon;
using namespace std;

typedef ListGraph::Node Node;
typedef ListGraph::Edge Edge;
typedef ListGraph::EdgeMap<int> LengthMap;
typedef lemon::ListGraph Graph;

#include <string>
#include <map>
#include <limits>

#include <iostream>
#include <fstream>
#include <lemon/lgf_reader.h>

void save_vector_int(std::vector<int> vec_to_save, std::string filename) {
    ofstream vector_file;
    vector_file.open(filename.c_str());
    //partition_file << partition_map.size() << endl;
    for (std::vector<int>::iterator itr = vec_to_save.begin(); itr != vec_to_save.end(); ++ itr) {
        vector_file << *itr << endl;
    }
}

void save_vector(std::vector<float> vec_to_save, std::string filename) {
    ofstream vector_file;
    vector_file.open(filename.c_str());
    //partition_file << partition_map.size() << endl;
    for (std::vector<float>::iterator itr = vec_to_save.begin(); itr != vec_to_save.end(); ++ itr) {
        vector_file << *itr << endl;
    }
}

void save_stabilities(umap &partitions_map, std::string filename) {
    ofstream partition_file;
    partition_file.open(filename.c_str());
    //partition_file << partition_map.size() << endl;
    for (umap::const_iterator part_itr = partitions_map.begin(); part_itr != partitions_map.end(); ++part_itr ) {
        partition_file << part_itr->second->key << ",";
        for (std::vector<double>::iterator st_itr = part_itr->second->stability.begin(); st_itr != part_itr->second->stability.end(); ++st_itr) {
            partition_file << *st_itr << ",";
        }
        partition_file << std::endl;
    }
    partition_file.close();
}


void load_expmats(std::string filename, std::vector<std::vector<double> > &expmats) {
    ifstream expmat_file(filename.c_str());
    std::string line;
    std::string expij;

    if ( !expmat_file.is_open() ) {
        std::cout << "couldn't open file" << endl;
        exit(1);
    }
    while (std::getline(expmat_file,line)) {
        std::stringstream lineStream(line);
        std::vector<double> current_expmat;
        while (std::getline(lineStream,expij,',')) {
            current_expmat.push_back(atof(expij.c_str()));
        }
        expmats.push_back(current_expmat);
    }
    expmat_file.close();
}

//USAGE ./find_parts filename(looks in data folder)
//Output: basname.stabilities, basename.taus, basename_landscape.lgf
int main(int argc, char* argv[])
{
    std::string UNIPROJ_HOME = std::getenv ("UNIPROJ_HOME");

    if (argc > 2) {
        cout << "Using Stability in range: t=" << atof(argv[2]) << " t < " << atof(argv[3]) << " t = t + t/" << atof(argv[4]) << std::endl;
        for (float t = atof(argv[2]); t<atof(argv[3]); t=t+t/atof(argv[4])) {
            m_times.push_back(t);
            //cout << t << endl;
        }
    }
    else {
        cout << "Using Modularity" << endl;
        m_times.push_back(1.0);
    }
    partitions_file.open("out.dat");
    partitions_vec.reserve(3000000);
    /*   graph *g; */
    unsigned long comp_tmp,neighbours,forbidden;
    int i;

    FILE * pFile;
    std::string basename = argv[1];
    load_expmats(UNIPROJ_HOME+"/data/expmats/"+basename+".expmat", expmats);

    //Find m and n
    /*for (std::vector<std::vector<double> >::iterator titr = expmats.begin(); titr != expmats.end(); ++titr) {
        double m_temp
        for (std::vector<double>::iterator eitr = itr->begin(); eitr != itr->end(); ++ eitr) {
             m_temp += *eitr;
        }
    }*/

    string filename = UNIPROJ_HOME+"/data/graphs/"+basename+".nke";
    cout << "Loading " << filename << endl;
    pFile = fopen (filename.c_str(),"r");

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
    /*  printf("begin\n"); */
    gen_partitions_prune_aux(comp_tmp, neighbours, forbidden);
    partitions_file.close();
    printf("%d\n",num_partitions);
    float avg_num_sets = num_sets / float(num_partitions);
    printf("%f\n",avg_num_sets);



    ListGraph g_landscape;
    for (i = 0;i<num_partitions;++i) {
        g_landscape.addNode();
    }
    int num_iterations = 0;

    cout << "saving markov times" << endl;
    save_vector(m_times,UNIPROJ_HOME+"/data/taus/"+basename+".taus");
    cout << "saving stabilities" << endl;
    save_stabilities(partition_map,UNIPROJ_HOME+"/data/stabilities/"+basename+"_unsorted.stabilities");

    //Find best (number of) partitions for each markov time
    std::vector<int> all_best_num_parts;
    std::vector<float> all_best_stabilities;
    for (unsigned int i =0;i<m_times.size(); ++i){
        double best_stability = std::numeric_limits<int>::min();
        int best_num_parts = 0;
        for (umap::const_iterator part_itr = partition_map.begin(); part_itr != partition_map.end(); ++part_itr ) {
            if (part_itr->second->stability[i] > best_stability) {
                best_stability = part_itr->second->stability[i];
                std::map<int,int> current_partition = bitset_to_map(part_itr->second->bit_partitions);
                best_num_parts = (std::max_element(current_partition.begin(), current_partition.end(),value_comparer))->second + 1;
            }
        }
        all_best_stabilities.push_back(best_stability);
        all_best_num_parts.push_back(best_num_parts);
    }
    save_vector(all_best_stabilities, UNIPROJ_HOME+"/data/stabilities/"+basename+"_best.stabilities");
    save_vector_int(all_best_num_parts, UNIPROJ_HOME+"/data/stabilities/"+basename+".best_nums");

    //return 0;

    //Find Neighbours
    for (umap::const_iterator part_itr = partition_map.begin(); part_itr != partition_map.end(); ++part_itr ) {

        if (num_iterations % 1000 == 0) {
            std::cout << "the number of iterations is" << num_iterations << std::endl;
        }
        uset partition_neighbours;
        find_neighbours(part_itr, partition_map, partition_neighbours);
        Node node_u = g_landscape.nodeFromId(part_itr->second->key);
        for (uset::iterator itr = partition_neighbours.begin(); itr != partition_neighbours.end();++itr) {
            Node node_v = g_landscape.nodeFromId((*itr)->second->key);
            Edge e = findEdge(g_landscape,node_u,node_v);
            if (e == INVALID) {
                Edge e = g_landscape.addEdge(node_u,node_v);
            }
        }
        num_iterations++;
    }
    cout << "Number of edges is " << countEdges(g_landscape) << " and is_connected is " << connected(g_landscape) << endl;
    graphWriter(g_landscape, UNIPROJ_HOME+"/data/landscapes/"+basename+"_landscape.lgf").run();

    free_graph(g);
    return 0;
}
