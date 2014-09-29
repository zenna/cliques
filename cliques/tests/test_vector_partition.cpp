#include <cliques/structures/vector_partition.h>
#include <cliques/tests/catch.hpp>

TEST_CASE( "VectorPartitions behaves correctly", "[VectorPartition]" ) {
  clq::VectorPartition A = clq::unassigned_partition(3);
  REQUIRE(A.node_count() == 3);
  REQUIRE(A.set_count() == 0);

  A = clq::global_partition(3);
  REQUIRE(A.set_count() == 1);  
  REQUIRE(A.find_set(0) == 0);
  REQUIRE(A.find_set(1) == 0);
  REQUIRE(A.find_set(2) == 0);

  A = clq::singletons_partition(3);
  REQUIRE(A.node_count() == 3);
  REQUIRE(A.set_count() == 3);
  REQUIRE(A.find_set(0) == 0);
  REQUIRE(A.find_set(1) == 1);
  REQUIRE(A.find_set(2) == 2);

  clq::unassign_node(A,2);
  REQUIRE(A.set_count() == 2);
  A.add_node_to_set(2,10);
  REQUIRE(A.set_count() == 3);
  REQUIRE(A.find_set(2) == 10);

  std::vector<int> v = {2,1,2,6};
  std::vector<int> v_norm = {0,1,0,2};
  clq::normalise_ids(v,false);
  REQUIRE(v == v_norm);

  clq::VectorPartition B{3, {1,4,4}, false};
  clq::VectorPartition C{3, {4,3,3}, false};
  clq::VectorPartition D{3, {4,3,4}, false};
  REQUIRE(B == C);
  REQUIRE(B != D);
}