#ifndef CLIQUES_NEIGHBOURS_H
#define CLIQUES_NEIGHBOURS_H

namespace cliques
{
/**
@brief  Louvain method - greedy algorithm to find community structure of a network.
@param[in]  my_graph     graph to find partition of
@param[in]  quality_function     partition quality function object
*/
template <typename P>
std::set<P> find_neighbours_group_moveset(P &partition)
{
	ads
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

}

#endif
