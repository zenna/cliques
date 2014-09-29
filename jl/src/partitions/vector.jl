#  @brief  A partition where the position in the vector denotes the node_id
#  and its value determines its set id.
#  E.g. [0,1,0] denotes the 0th and 2nd node belong to the same set (0),
#  and the 1st node belongs to its own set.
#  Good when you need constant time assignment of a node to a set
#  and removal of a node from a set

type VectorPartition
  n_nodes::Int64
  set_ids::Vector{Int64}
  is_normalised::Bool
end

copy(p::VectorPartition) = VectorPartition(p.n_nodes,p.set_ids,p.is_normalised)

const no_set_id = -1 #TODO This really doesn't make sense in the context of a partition, can we get rid of it?

## ===========
## Constructors
singleton_partition(n_nodes::Int64) = VectorPartition(n_nodes, collect(0:n_nodes-1), true)
global_partition(n_nodes::Int64) = VectorPartition(n_nodes, zeros(Int64, n_nodes), true)
unassigned_partition(n_nodes::Int64) = VectorPartition(n_nodes, fill(no_set_id, n_nodes), false)

find_set(p::VectorPartition, node_id::Int64) = p.set_ids[node_id]
add_node_to_set!(p::VectorPartition, node_id, set_id) = p.set_ids[node_id] = set_id
unassign_node!(p::VectorPartition, node_id::Int64) = p.set_ids[node_id] = no_set_id
function normalise_ids!(p::VectorPartition)
  set_new_to_old = Dict{Int64,Int64}()
  set_old_to_new = Dict{Int64,Int64}()
  if !p.is_normalised
    new_setid = 0
    for i = 1:length(p.set_ids)
      set_id = p.set_ids[i]
      if haskey(set_old_to_new, set_id)
        p.set_ids[i] = set_old_to_new[set_id]
      else
        set_new_to_old[new_setid] = set_id
        set_old_to_new[set_id] = new_setid
        p.set_ids[i] = new_setid
        new_setid += 1
      end
    end
  else
    for set_id in p.set_ids
      set_new_to_old[set_id] = set_id
    end
  end
  p.is_normalised = true
  set_new_to_old
end

node_count(p::VectorPartition) = length(p.set_ids)
function unique_sets(p::VectorPartition)
  seen = Set()
  for set_id in p.set_ids
    if set_id != no_set_id
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
  for i = 1:length(p.set_ids)
    if p.set_ids[i] == set_id
      push!(nodes,i)
    end
  end
end

function ==(p::VectorPartition, q::VectorPartition)
  if p.is_normalised && q.is_normalised
    p.set_ids == q.set_ids
  else
    pcopy = copy(p)
    qcopy = copy(q)
    normalise_ids!(pcopy)
    normalise_ids!(qcopy)
    p.set_ids == q.set_ids
  end
end

is_assigned(p::VectorPartition, node_id::Int64) = p.set_ids[node_id] != no_set_id
is_unassigned(p::VectorPartition, node_id::Int64) = p.set_ids[node_id] == no_set_id
