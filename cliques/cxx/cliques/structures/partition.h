#ifndef CLIQUES_PARTITION_H
#define CLIQUES_PARTITION_H

#include <bitset>
#include <set>
#include <cliques/structures/graph.h>

namespace cliques {
    enum new_partition_type {
        GLOBULAR,
        SINGLETONS
    };

    struct Partition // the event log data structure
    {
        int key;
        std::bitset<128> bit_partitions;
        double modularity;
        std::vector<double> stability;
    };
} //namespace cliques

#endif //CLIQUES_STRUCTURES_H
