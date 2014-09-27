#  @brief  A partition where the position in the vector denotes the node_id
#  and its value determines its set id.
#  E.g. [0,1,0] denotes the 0th and 2nd node belong to the same set (0),
#  and the 1st node belongs to its own set.
#  Good when you need constant time assignment of a node to a set
#  and removal of a node from a set

type VectorPartition
  n_nodes::Int64
  partition_vector::Vector{Int64}
  is_normalised::Bool
end

const nosetid = -1 #TODO This really doesn't make sense in the context of a partition, can we get rid of it?
collect(0:3)
## ===========
## Constructors
singletons(n_nodes::Int64) = VectorPartition(n_nodes, collect(0:n_nodes-1), true)
globalpartition(n_nodes::Int64) = VectorPartition(n_nodes, zeros(Int64, n_nodes), true)
unassignedpartition(n_nodes::Int64) = VectorPartition(n_nodes, fill(nosetid, n_nodes), false)

find_set(p::VectorPartition, node_id::Int64) = p.partition_vector[node_id]
add_node_to_set!(p::VectorPartition, node_id, set_id) = p.partition_vector[node_id] = set_id
unassign_node!(p::VectorPartition, node_id::Int64) = p.partition_vector[node_id] = nosetid
function normalise_ids!(p::VectorPartition)
  set_new_to_old = Dict{Int64,Int64}()
  set_old_to_new = Dict{Int64,Int64}()
  if !p.is_normalised
    new_setid = 0
    for i = 1:length(p.partition_vector)
      set_id = p.partition_vector[i]
      if haskey(set_old_to_new, set_id)
        p.partition_vector[i] = set_old_to_new[set_id]
      else
        set_new_to_old[new_setid] = set_id
        set_old_to_new[set_id] = new_setid
        p.partition_vector[i] = new_setid
        new_setid += 1
      end
    end
  else
    for set_id in p.partition_vector
      set_new_to_old[set_id] = set_id
    end
  end
  p.is_normalised = true
  set_new_to_old
end

node_count(p::VectorPartition) = length(p.partition_vector)
function unique_sets(p::VectorPartition)
  seen = Set()
  for set_id in p.partition_vector
    if set_id != nosetid
      push!(seen, set_id)
    end
  end
  seen
end

function set_count(p::VectorPartition)
  length(unique_sets(p))
end

function nodes_from_set(p::VectorPartition, set_id::Integer)
  nodes = Int64[]
  for i = 1:length(p.partition_vector)
    if p.partition_vector[i] == set_id
      push!(nodes,i)
    end
  end
end

function ==(p::VectorPartition, q::VectorPartition)
  if p.is_normalised && q.is_normalised
    p.partition_vector == q.partition_vector
  else
    pcopy = copy(p)
    qcopy = copy(q)println()
    normalise_ids!(pcopy)
    normalise_ids(qcopy)
    p.partition_vector == q.partition_vector
  end
end

is_assigned(p::VectorPartition, node_id::Int64) = p.partition_vector[node_id] != nosetid
is_unassigned(p::VectorPartition, node_id::Int64) = p.partition_vector[node_id] == nosetid
