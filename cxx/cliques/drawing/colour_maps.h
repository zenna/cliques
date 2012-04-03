#ifndef CLIQUES_COLOURMAPS_H
#define CLIQUES_COLOURMAPS_H

#include <map>
#include <vector>

#include <cliques/algorithms/scalers.h>
#include <cliques/structures/common.h>

namespace cliques {

/**
 @brief  This colour map functor, gives nodes colours based on their group
 */
template<typename P>
struct make_partition_colour_map {
    P &partition_;

    make_partition_colour_map(P &partition) :
        partition_(partition) {
    }

    std::map<int, cliques::xyz_colour<float> > operator ()() {
        std::map<int, cliques::xyz_colour<float> > colour_map;
        int num_sets = partition_.set_count();
        float hue = 0.0;
        float hue_increment = (1.0 - hue) / num_sets;
        for (typename P::PartIterator pitr = partition_.begin(); pitr
                != partition_.end(); ++pitr) {
            for (typename P::NodeIterator nitr = pitr.begin(); nitr
                    != pitr.end(); ++nitr) {
                // insert new colour into map
                colour_map.insert(std::pair<int, cliques::xyz_colour<float> >(
                        *nitr, cliques::xyz_colour<float>(hue, 0.7, 0.7)));
            }
            hue = hue + hue_increment;
        }
        return colour_map;
    }
};

/**
 @brief  This colour map functor, gives nodes colours based on their group
 */
template<typename P>
struct highlight_group_node_map {
    P &partition_;
    int highlighted_set;

    highlight_group_node_map(P &partition, int highlighted_set) :
        partition_(partition), highlighted_set(partition_.find_set(
                highlighted_set)) {
    }

    std::map<int, cliques::xyz_colour<float> > operator ()() {
        std::map<int, cliques::xyz_colour<float> > colour_map;

        // Inefficient!
        for (typename P::PartIterator pitr = partition_.begin(); pitr
                != partition_.end(); ++pitr) {
            for (typename P::NodeIterator nitr = pitr.begin(); nitr
                    != pitr.end(); ++nitr) {
                // insert new colour into map
                if (partition_.find_set(*nitr) == highlighted_set) {
                    colour_map.insert(
                            std::pair<int, cliques::xyz_colour<float> >(*nitr,
                                    cliques::xyz_colour<float>(1.0, 0.7, 0.7)));
                } else {
                    colour_map.insert(
                            std::pair<int, cliques::xyz_colour<float> >(*nitr,
                                    cliques::xyz_colour<float>(0.0, 0.0, 0.3)));
                }
            }
        }
        return colour_map;
    }
};

/**
 @brief  This colour map functor gives nodes colours based on their partition group
 */
struct make_energy_edge_colour_map {
    std::vector<float> &energies_;
    float max_weight;
    float min_weight;
    bool DYNAMIC_SCALE;

    make_energy_edge_colour_map(std::vector<float> &energies) :
        energies_(energies), max_weight(-100000.00), min_weight(100000.00),
                DYNAMIC_SCALE(true) {
    }

    std::map<int, cliques::xyz_colour<float> > operator ()() {
        std::map<int, cliques::xyz_colour<float> > colour_map;
        //assuming there is a mapping between node_id and position in energies vector

        // Find range
        if (DYNAMIC_SCALE == true) {
            for (std::vector<float>::iterator itr = energies_.begin(); itr
                    != energies_.end(); ++itr) {
                if (*itr > max_weight) {
                    max_weight = *itr;
                }
                if (*itr < min_weight) {
                    min_weight = *itr;
                }
            }
        }

        float shift = 0.0 - min_weight;
        float range = max_weight - min_weight;
        float hue;

        // TODO account for case when range == 0
        // currently causes nan due to divide by 0
        int i = 0;
        for (std::vector<float>::iterator itr = energies_.begin(); itr
                != energies_.end(); ++itr) {

            hue = (*itr + shift) / range;

            colour_map.insert(std::pair<int, cliques::xyz_colour<float> >(i,
                    cliques::xyz_colour<float>(hue, 0.7, 0.7)));
            ++i;
        }
        return colour_map;
    }
};

}

#endif
