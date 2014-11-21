using Cliques
import Cliques: unassignedpartition, globalpartition, set_count, find_set, node_count, add_node_to_set!,
                normalise_ids!, singletons, unassign_node!
using Base.Test

A = unassigned_partition(3)
@test node_count(A) == 3
@test set_count(A) == 0

A = global_partition(3)
@test set_count(A) == 1
@test find_set(A,1) == 0
@test find_set(A,2) == 0
@test find_set(A,3) == 0

A = singletons_partition(3)
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
@test A.set_ids == [0,1,0,2]

A = VectorPartition(3, [1,4,4], false)
B = VectorPartition(3, [4,3,3], false)
C = VectorPartition(3, [4,3,4], false)
@test A == B
@test A != C
