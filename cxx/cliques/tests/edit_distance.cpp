#include <iostream>
#include <cliques/structures/vector_partition.h>
#include <cliques/helpers/edit_distance.h>
#include <cliques/helpers.h>

int main() {
    cliques::VectorPartition A(5);
    A.add_node_to_set(0,0);
    A.add_node_to_set(1,0);
    A.add_node_to_set(2,1);
    A.add_node_to_set(3,1);
    A.add_node_to_set(4,2);
    cliques::VectorPartition B(5);
    B.add_node_to_set(0,0);
    B.add_node_to_set(1,1);
    B.add_node_to_set(2,2);
    B.add_node_to_set(3,2);
    B.add_node_to_set(4,2);

    cliques::Hungarian hungarian(A,B);
    int edit_distance = hungarian.edit_distance();
    cliques::output(edit_distance);
    return 0;
};
