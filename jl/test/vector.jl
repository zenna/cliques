using Cliques
import Cliques: unassignedpartition, globalpartition, set_count, find_set, node_count, add_node_to_set!,
                normalise_ids!, singletons, unassign_node!
using Base.Test

A = unassignedpartition(3)
@test node_count(A) == 3
@test set_count(A) == 0

A = globalpartition(3)
@test set_count(A) == 1
@test find_set(A,1) == 0
@test find_set(A,2) == 0
@test find_set(A,3) == 0

A = singletons(3)
@test node_count(A) == 3
@test set_count(A) == 3
@test find_set(A,1) == 0
@test find_set(A,2) == 1
@test find_set(A,3) == 2

unassign_node!(A,2)
@test set_count(A) == 2
add_node_to_set!(A,2,10)
@test set_count(A) == 3
@test find_set(A,2) == 10

A = VectorPartition(4,[2,1,2,6], false)
normalise_ids!(A)
@test A.partition_vector == [0,1,0,2]