/* Copyright (c) Z Tavares, M Schaub - 2010-2014 */
#pragma once

#include <vector>
#include <set>
#include <map>

namespace clq {

// A node may not be yet assigned a set id
const int no_set_id = -1;

/**
 @brief  A partition where the position in the vector denotes the node_id
 and its value determines its set id.

 E.g. [0,1,0] denotes the 0th and 2nd node belong to the same set (0),
 and the 1st node belongs to its own set.

 Good when you need constant time assignment of a node to a set
 and removal of a node from a set
 */
class VectorPartition {
private:
  int num_nodes;
  std::vector<int> set_ids;
  bool is_normalised_;

public:
  //#################### CONSTRUCTORS ####################

  explicit VectorPartition(int num_nodes, std::vector<int> set_ids,
                           bool is_normalised_) :
    num_nodes(num_nodes),
    set_ids(set_ids),
    is_normalised_(is_normalised_) {
  }

  // construct partition from vector
  explicit VectorPartition(std::vector<int> partition) :
    num_nodes(partition.size()), set_ids(partition), is_normalised_(false) {
  }

  //#################### PUBLIC METHODS ####################
public:
  /** @brief  set_id assigned to node  */
  int find_set(int node_id) const {
    return set_ids[node_id];
  }

  // Assign a node to a set
  void add_node_to_set(int node_id, int set_id) {
    set_ids[node_id] = set_id;
  }


  /** @brief  Number of elements in partition  */
  int node_count() const { return set_ids.size(); }

  /** @brief  Number of subsets (groups) in partition  */
  int set_count() const {
    std::set<int> seen_nodes;
    for (int set_id : set_ids) {
      if (set_id != -1) {
        seen_nodes.insert(set_id);
      }
    }
    return seen_nodes.size();
  }

  friend std::map<int, int> normalise_ids(VectorPartition &p);
  
  // whether partition is normalised is stored for efficiency.
  bool is_normalised() const {return is_normalised_;}
  void normalise() {is_normalised_ = true;}
  void unnormalise() {is_normalised_ = false;}

  friend bool operator==(const VectorPartition& lhs,
                         const VectorPartition& rhs);
  friend bool operator!=(const VectorPartition& lhs,
                         const VectorPartition& rhs);
};

// =====================
// External Constructors

/** @brief  Partition where every node assigned to its own group  */
VectorPartition singletons_partition(int n_nodes) {
  unsigned int i = 0;
  std::vector<int> set_ids(n_nodes);
  for (int &x : set_ids) {
    x = i;
    i += 1;
  }
  return VectorPartition(n_nodes, set_ids, true);
}

/** @brief  Partition with all elements in same group (0)  */
VectorPartition global_partition(int n_nodes) {
  std::vector<int> set_ids(n_nodes, 0);
  return VectorPartition(n_nodes, set_ids, true);
}

/** @brief  Partition with no nodes assigned to a group  */
VectorPartition unassigned_partition(int n_nodes) {
  std::vector<int> set_ids(n_nodes, no_set_id);
  return VectorPartition(n_nodes, set_ids, false);
}

/** @brief  all nodes with given set_id  */
std::vector<int> nodes_from_set(VectorPartition p, int set_id) {
  std::vector<int> nodes_in_set;
  for (int i = 0; i < p.node_count(); ++i) {
    if (p.find_set(i) == set_id) {
      nodes_in_set.push_back(i);
    }
  }
  return nodes_in_set;
}

/** @brief  remove node from partition: remove group assignment */
void unassign_node(VectorPartition &p, int node_id) {
  p.add_node_to_set(node_id, no_set_id);
  p.unnormalise();
}

/** @brief  is given node_id assigned a group?  */
bool const is_assigned(VectorPartition p, int node_id) {
  return p.find_set(node_id) != no_set_id;
}

/** @brief  is given node_id not assigned a group?  */
bool const is_unassigned(VectorPartition p, int node_id) {
  return p.find_set(node_id) == no_set_id;
}

std::map<int, int> normalise_ids(std::vector<int> &p, bool is_normalised) {
  std::map<int, int> set_new_to_old;
  std::map<int, int> set_old_to_new;

  if (!is_normalised) {
    int new_set_id = 0;
    int i = 0;
    for (int set_id : p) { //FIXME: NO ACCESS TO PRIMVATE
      auto old_set = set_old_to_new.find(set_id);
      if (old_set != set_old_to_new.end()) {
        p[i] = old_set->second;
      }
      else {
        set_new_to_old[new_set_id] = set_id;
        set_old_to_new[set_id] = new_set_id;
        p[i] = new_set_id;
        new_set_id += 1;
      }
      i += 1;
    }
  }
  else {
    for (int set_id : p) {
      set_new_to_old[set_id] = set_id;
    }
  }
  return set_new_to_old;
}

/** @brief normalises a partition to a canonical form

    normalised vectors assign set_ids contiguously starting at 0
    with 0th node, and increasing only with increasing node_id */
std::map<int, int> normalise_ids(VectorPartition &p) {
  auto set_new_to_old = normalise_ids(p.set_ids, p.is_normalised());
  p.normalise();
  return set_new_to_old;
}

/** @brief Semantic equivalence, are the groupings identical? */
bool operator==(const VectorPartition& lhs, const VectorPartition& rhs) {
  if (lhs.is_normalised() && rhs.is_normalised()) {
    return (lhs.set_ids == rhs.set_ids);
  }
  else {
    VectorPartition a(lhs);
    VectorPartition b(rhs);

    normalise_ids(a);
    normalise_ids(b);
    return a.set_ids == b.set_ids;
  }
}

bool operator!=(const VectorPartition& lhs, const VectorPartition& rhs) {
  return !(lhs == rhs);
}

}