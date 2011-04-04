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
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>

#include <lemon/list_graph.h>
#include <lemon/connectivity.h>
#include <lemon/concepts/graph_components.h>
#include <lemon/concepts/graph.h>

//#include <boost/functional/hash.hpp>
//#include <boost/unordered_map.hpp>
//#include <boost/unordered_set.hpp>

#define MAX_LINE_LENGTH 1000
#define SET_MAX_SIZE 63

#define N 128
#define NUM_NODES 16
#define BITS_PER_NODE 4
using namespace std;

//Global Vars
typedef std::vector<std::vector<int> > vec_2d_ints;
std::map<std::pair<int,int>,float> w_weights;
std::map<int,float> w_degrees;
float w_m;

void reorder_ints(int &i, int &j) {
    int buffer;
    if (i <= j) {
        return;
    }
    else {
        buffer = i;
        i = j;
        j = buffer;
    }

}

int reorder_map(std::map<int,int> &partition_map) {
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
    return start_part;
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


typedef struct graph
{
    int n;
    int m;
    int **links;
    int *degrees;
    int *capacities;
} graph;

/* Global variables */
graph *g;

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

float A(int one, int two) {
    for(int i=0; i<g->degrees[one]; i++){
        if (g->links[one][i] == two) {
            return 1.0;
        }
    }

    return 0.0;
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
    std::vector<std::vector<int > > p = map_to_vector(partition_map);

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

void find_stability(std::map<int,int> &partition_map, std::map<float,std::vector<float> > &stabilities, std::vector<float> markov_times) {
    std::vector<std::vector<int > > p = map_to_vector(partition_map);
    float two_m = 2*g->m;
    float first_term = 0.0;
    float second_term = 0.0;

    for (int s = 0;  s < p.size(); ++s) {
        for (int i = 0; i< p[s].size(); ++i) {
            for (int j = 0; j < p[s].size(); ++j) {
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
        stabilities[*t].push_back(R);
    }
}

float W(int one, int two) {
    //Look up these weights in a
    reorder_ints(one,two);
    return w_weights[make_pair(one,two)];
}

void load_weights(std::string filename, std::map<std::pair<int,int>,float> &w_weights, std::map<int,float> &w_degrees, float &w_m) {
    w_m = 0.0;
    //Create Pair
    //std::map<std::pair<int,int>,float> w_weights;

    //Degrees
    // Go through each line counting
        // if i
    ifstream weights_file(filename.c_str());
    std::string line;
    std::string weight;

    if ( !weights_file.is_open() ) {
        std::cout << "couldn't open file" << endl;
        exit(1);
    }
    int counter = 0;
    int sample_in_array = 0;
    int i = 0;
    while (std::getline(weights_file,line)) {
        std::stringstream lineStream(line);

        int j = 0;
        while (std::getline(lineStream,weight,' ')) {
            if (i != j) {
                w_degrees[i] += atof(weight.c_str());
                w_m += atof(weight.c_str());
                int t_i = i;
                int t_j = j;
                reorder_ints(t_i,t_j);
                std::pair<int,int>edge (t_i,t_j);
                w_weights[edge] = atof(weight.c_str());
            }
            j++;
        }
        i++;
    }
    weights_file.close();
}

void find_stability_weighted(std::map<int,int> &partition_map, std::map<float,std::vector<float> > &stabilities, std::vector<float> markov_times) {
    std::vector<std::vector<int > > p = map_to_vector(partition_map);
    float two_m = 2*w_m;
    float first_term = 0.0;
    float second_term = 0.0;

    for (int s = 0;  s < p.size(); ++s) {
        for (int i = 0; i< p[s].size(); ++i) {
            for (int j = 0; j < p[s].size(); ++j) {
                first_term = first_term + float(w_degrees[p[s][i]]*w_degrees[p[s][j]]);
                second_term = second_term + W(p[s][i],p[s][j]);
            }
        }
    }

    first_term = first_term/(two_m*two_m);
    second_term = second_term/two_m;

    float R;
    for (std::vector<float>::iterator t = markov_times.begin(); t != markov_times.end(); ++t) {
        R = 0.0;
        R = (1-*t) - first_term + *t * second_term;
        stabilities[*t].push_back(R);
    }
}

std::vector<float> m_times;

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

using namespace lemon;
using namespace std;

typedef ListGraph::Node Node;
typedef ListGraph::Edge Edge;
typedef ListGraph::EdgeMap<int> LengthMap;
typedef lemon::ListGraph Graph;

#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <lemon/lgf_reader.h>
#include <lemon/random.h>


void save_vector(std::vector<float> vec_to_save, std::string filename) {
    ofstream vector_file;
    vector_file.open(filename.c_str());
    //partition_file << partition_map.size() << endl;
    for (std::vector<float>::iterator itr = vec_to_save.begin(); itr != vec_to_save.end(); ++ itr) {
        vector_file << *itr << endl;
    }
}

std::map<float,std::vector<float> > random_walk(int num_steps, bool find_weighted_stability = false) {
    std::map<int,int> current_partition;
    //int num_parts = rnd((g->n) - 1);
    int num_parts;
    bool is_first_map_connected = false;
    while (is_first_map_connected == false) {
        for (int i = 0;i<g->n;++i) {
            current_partition[i] = rnd.integer(g->n);
        }
        num_parts = reorder_map(current_partition);
        vec_2d_ints partition_vec = map_to_vector(current_partition);
        is_first_map_connected = true;
        for (int i =0;i<partition_vec.size();++i) {
            if (is_connected(partition_vec[i]) == false) {
                is_first_map_connected = false;
                break;
            }
        }
        //print_map(current_partition);
    }

    //Check connected
    //find_stability(partition_map, stabilities, m_times);
    std::map<float,std::vector<float> > stabilities;

    for (int i=0;i<num_steps;++i) {
        //Find random neighbour
        std::map<int,int> neighbour_partition;
        //print_map(current_partition);
        bool is_map_connected = false;
        int node_to_move;
        do {
            neighbour_partition = current_partition;
            //node_to_move;
            int current_part;
            int next_part;

            do {
                node_to_move = rnd.integer(g->n);
                current_part = neighbour_partition[node_to_move];
                next_part = rnd.integer(num_parts+1);
            }
            while (current_part == next_part);
            neighbour_partition[node_to_move] = next_part;


            //Check Connected

            // More efficient than checking whole graph (I think?)
            std::vector<int> first_changed_part, last_changed_part;
            for (std::map<int,int>::iterator new_itr = neighbour_partition.begin(); new_itr != neighbour_partition.end(); ++new_itr) {
                if (new_itr->second == current_part)
                    first_changed_part.push_back(new_itr->first);
                if (new_itr->second == next_part)
                    last_changed_part.push_back(new_itr->first);
            }
            is_map_connected = ((!first_changed_part.size()) || is_connected(first_changed_part) ) && ( (!last_changed_part.size()) || is_connected(last_changed_part) );
        }
        while (is_map_connected == false);
        //cout << "unordered neighbour is ";
        //print_map(neighbour_partition);
        num_parts = reorder_map(neighbour_partition);
        //cout << "current is ";
        //print_map(current_partition);
        //cout << "reordered neighbour is ";
        //print_map(neighbour_partition);
        if (neighbour_partition == current_partition) {
            //print_map(neighbour_partition);
            i--;
            continue;
        }
        current_partition = neighbour_partition;
        //print_map(current_partition);
        if (find_weighted_stability == true) {
            find_stability_weighted(current_partition, stabilities, m_times);
        }
        else {
            find_stability(current_partition, stabilities, m_times);
        }
    }

    return stabilities;
}

std::vector<int> find_symbols(std::vector<float> &stabilities, float eta) {
    std::vector<int> symbols;
    for (std::vector<float>::iterator old_stab = stabilities.begin(); old_stab != (stabilities.end()-1); ++old_stab) {
        float stab_diff = *(old_stab+1) - *old_stab;
        if (stab_diff < -eta) {
            symbols.push_back(-1);
        }
        else if (abs(stab_diff) <= eta) {
            symbols.push_back(0);
        }
        else if (stab_diff  > eta) {
            symbols.push_back(1);
        }
    }
    return symbols;
}

#include <math.h>
float calculate_entropy(std::vector<int> &symbols) {
    float H = 0.0;

    for (int p=-1;p<2;++p) {
        for (int q=-1;q<2;++q) {
            if (p!=q) {
                float P = 0;
                for (unsigned int i=0;i<symbols.size()-1 ;++i) {
                    if ((symbols[i] == p) && (symbols[i+1] == q)) {
                        P++;
                    }
                }
                P = P/symbols.size();
                float logP = log(P)/log(6);
                if (P!=0)
                    H = H + P*logP;
            }
        }
    }
    return -H;
}

//USAGE ,/random_walk graph_filename weight_filename start_stability_index end_stability_index
int main(int argc, char* argv[]) {
    std::string UNIPROJ_HOME = std::getenv ("UNIPROJ_HOME");

    bool do_weighted_analysis = false;
    if (argc == 3) {
        //cout << "doing weighted analysis" << endl;
        string filename = UNIPROJ_HOME+"/data/weights/";
        filename += argv[2];
        load_weights(filename,w_weights, w_degrees, w_m);
        do_weighted_analysis = true;
    }

    for (float t = 0.008; t<10000.0; t=t+t/25) {
        m_times.push_back(t);
    }
    FILE * pFile;

    //pFile = fopen(graph.c_str(), "r");
    std::string basename = argv[1];

    string filename = UNIPROJ_HOME+"/data/graphs/"+basename+".nke";
    pFile = fopen (filename.c_str(),"r");
    g=graph_from_file(pFile);

    //int num_walks = 30;
    int num_steps = 30000;
    std::map<float,std::vector<float> > alpha = random_walk(num_steps, do_weighted_analysis);

    // Calculate ETas
    std::map<float,float> eta_stars;
    for (std::map<float,std::vector<float> >::iterator itr = alpha.begin(); itr != alpha.end(); ++itr) {
        eta_stars[itr->first] = 0.0;
        for (std::vector<float>::iterator fitr = itr->second.begin(); fitr != (itr->second.end()-1); ++fitr) {
            float step = abs(*fitr - *(fitr+1));//<< endl;
            if (step > eta_stars[itr->first]) {
                eta_stars[itr->first] = step;
            }
        }
    }

    // Calculate Entropies
    int num_eta_steps = 100;
    std::map<float,float> entropies;
    for (std::map<float,std::vector<float> >::iterator itr = alpha.begin(); itr != alpha.end(); ++itr) {
        entropies[itr->first] = 0.0;
        for (float eta = 0.0; eta < eta_stars[itr->first]; eta = eta + eta_stars[itr->first]/num_eta_steps) {
            std::vector<int> symbols = find_symbols(itr->second,eta);
            float entropy = calculate_entropy(symbols);
            if (entropy > entropies[itr->first]) {
                entropies[itr->first] = entropy;
            }
        }
    }

    for (std::map<float,float>::iterator itr = entropies.begin(); itr != entropies.end(); ++itr) {
        cout << itr->second << endl;
    }


    free_graph(g);
    return 0;
}

