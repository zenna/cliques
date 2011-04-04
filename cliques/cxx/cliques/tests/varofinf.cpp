#include <iostream>
#include <helpers.h>
#include <structures/disjointset.h>
#include <algorithms/varofinf.h>

int main() {
	cliques::DisjointSetForest<int> partition1;
    int num_elements = 10;
    for (int i =0; i<num_elements; ++i) {
    	partition1.add_element(i);
    }

    partition1.union_sets(0,1);
    partition1.union_sets(2,3);
    partition1.union_sets(3,4);
    partition1.union_sets(4,5);
    partition1.union_sets(9,8);

	cliques::DisjointSetForest<int> partition2;
    for (int i =0; i<num_elements; ++i) {
    	partition2.add_element(i);
    }

    partition2.union_sets(0,1);
    partition2.union_sets(6,7);
    partition2.union_sets(3,4);
    partition2.union_sets(1,3);
    partition2.union_sets(6,8);

    cliques::print_partition(partition1);
    std::cout << std::endl;
    cliques::print_partition(partition2);

    std::cout << "Variation of information is "
    		<< cliques::find_variation_of_information(partition1,partition2)
    		<< std::endl;

    return 0;
};
