#include <iostream>
#include <cliques/helpers.h>
#include <cliques/structures/disjointset.h>
#include <cliques/algorithms/varofinf.h>

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

    /*VecPartition test(6);
        VecPartition test2(6);

        std::vector<int> v = {1,2,2,3,0,1};
        std::vector<int> v2 = {5,3,3,2,1,5};

        test.partition_vector = v;
        test2.partition_vector = v2;

        bool same = test == test2;
        std::cout << "same" << same <<std::endl;

        */

    return 0;
};
