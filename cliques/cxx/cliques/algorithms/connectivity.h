#ifndef CLIQUES_NEIGHBOURS_H
#define CLIQUES_NEIGHBOURS_H

namespace cliques
{
/**
@brief  Louvain method - greedy algorithm to find community structure of a network.
@param[in]  my_graph     graph to find partition of
@param[in]  quality_function     partition quality function object
*/
template <typename P, typename G>
bool is_connected(P std::vector<int> part, G) {
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


}

#endif
