#include <iostream>
#include <helpers.h>
#include <structures/disjointset.h>

//TODO
//REMOVE NEIGHBOURS TO SELF
int main() {
	cliques::DisjointSetForest<int> partition;
    int num_elements = 10;
    for (int i =0; i<num_elements; ++i) {
    	partition.add_element(i);
    }

    partition.union_sets(0,1);
    partition.union_sets(2,3);
    partition.union_sets(3,4);
    partition.union_sets(4,5);
    partition.union_sets(9,8);

    cliques::print_partition(partition);

    for (cliques::DisjointSetForest<int>::PartIterator pitr = partition.begin();
    		pitr != partition.end(); ++pitr) {
    	std::cout << partition.set_size(pitr) << std::endl;
    }



    return 0;
};
