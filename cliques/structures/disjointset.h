/***
 * millipede: DisjointSetForest.h
 * Copyright Stuart Golodetz, 2009. All rights reserved.
 * Modifed by Zenna Tavares 2011
 ***/
#pragma once

#include <map>
#include <set>
#include <iostream>

#include <unordered_map>

namespace clq {

/**
 @brief  A disjoint set forest is a fairly standard data structure used to represent the partition of
 a set of elements into disjoint sets in such a way that common operations such as merging two
 sets together are computationally efficient.

 This implementation uses the well-known union-by-rank and path compression optimizations, which together
 yield an amortised complexity for key operations of O(a(n)), where a is the (extremely slow-growing)
 inverse of the Ackermann function.

 The implementation also allows clients to attach arbitrary data to each element, which can be useful for
 some algorithms.

 @tparam T   The type of data to attach to each element (arbitrary)
 */
template<typename T>
class DisjointSetForest {
	//#################### NESTED CLASSES ####################
public:
	struct Element {
		T m_value;
		int m_parent;
		int m_rank;

		Element(const T& value, int parent) :
			m_value(value), m_parent(parent), m_rank(0) {
		}
	};

	//#################### PRIVATE VARIABLES ####################
private:
	//mutable std::map<int,Element> m_elements;
	mutable std::unordered_map<int, Element> m_elements;
	unsigned int m_setCount;
	Element *last_unioned_node_;
	int last_unioned_parent_;
	bool was_last_union_effective;

	//#################### CONSTRUCTORS ####################
public:
	/**
	 @brief  Constructs an empty disjoint set forest.
	 */
	DisjointSetForest() :
		m_setCount(0) {
	}

	/**
	 @brief  Constructs a disjoint set forest from an initial set of elements and their associated values.

	 @param[in]  initialElements     A map from the initial elements to their associated values
	 */
	explicit DisjointSetForest(
			const std::unordered_map<int, T>& initialElements) :
		m_setCount(0) {
		add_elements(initialElements);
	}

	//#################### ITERATORS ####################
	class NodeIterator;

	class PartIterator {
	public:
		PartIterator(
				typename std::unordered_map<int, Element>::iterator part_itr,
				DisjointSetForest &parent_disjoint_set) :
			part_itr_(part_itr), parent_disjoint_set_(parent_disjoint_set) {
		}

		PartIterator& operator=(const PartIterator& other) {
			part_itr_ = other.part_itr;
			visited_parents_ = other.visited_parents_;
			return (*this);
		}

		bool operator==(const PartIterator& other) {
			return (part_itr_ == other.part_itr_);
		}

		bool operator!=(const PartIterator& other) {
			return (part_itr_ != other.part_itr_);
		}

		PartIterator& operator++() {
			int parent_set = parent_disjoint_set_.find_set(part_itr_->first);
			visited_parents_.insert(parent_set);

			if (visited_parents_.size() == parent_disjoint_set_.m_setCount) {
				part_itr_ = parent_disjoint_set_.m_elements.end();
			}

			else {
				do {
					++part_itr_;
					parent_set
							= parent_disjoint_set_.find_set(part_itr_->first);
				} while (visited_parents_.find(parent_set)
						!= visited_parents_.end());
				//visited_parents_.insert(parent_set);
			}

			return *this;
		}

		int operator*() {
			return (parent_disjoint_set_.find_set(part_itr_->first));
		}

		NodeIterator begin() {
			typename std::unordered_map<int, Element>::iterator
					tmp_node_itr_ = part_itr_;
			return NodeIterator(tmp_node_itr_, this->parent_disjoint_set_);
		}

		NodeIterator end() {
			typename std::unordered_map<int, Element>::iterator
					tmp_node_itr_ = parent_disjoint_set_.m_elements.end();
			return NodeIterator(tmp_node_itr_, this->parent_disjoint_set_);
		}

	private:
		typename std::unordered_map<int, Element>::iterator part_itr_;
		//typename std::map<int,T>::iterator part_itr_;
		std::set<int> visited_parents_;
		DisjointSetForest &parent_disjoint_set_;
		//Node* node_;
	};

	PartIterator begin() {
		//typename std::map<int,T>::iterator tmp_part_itr_;
		//tmp_part_itr_ = m_elements.begin();
		typename std::unordered_map<int, Element>::iterator tmp_part_itr_ =
				m_elements.begin();
		return PartIterator(tmp_part_itr_, *this);
	}

	PartIterator end() {
		typename std::unordered_map<int, Element>::iterator tmp_part_itr_ =
				m_elements.end();
		return PartIterator(tmp_part_itr_, *this);
	}

	class NodeIterator {
	public:
		NodeIterator(
				typename std::unordered_map<int, Element>::iterator node_itr,
				DisjointSetForest &parent_disjoint_set) :
			node_itr_(node_itr), parent_disjoint_set_(parent_disjoint_set) {
			if (node_itr_ != parent_disjoint_set.m_elements.end()) {
				parent_set_ = parent_disjoint_set_.find_set(node_itr_->first);
				//num_nodes_in_set_ = parent_disjoint_set_.element_count(parent_set_);
				//num_nodes_seen_ = 1;
			}
		}

		NodeIterator& operator++() {
			/*if (num_nodes_seen_ == num_nodes_in_set_) {
			 node_itr_ = parent_disjoint_set_.m_elements.end();
			 }*/

			// Iterate node_itr_ until we find a node of the same set
			// or we reach the end of all elements
			//int new_parent_set;
			do {
				++node_itr_;
			} while (node_itr_ != parent_disjoint_set_.m_elements.end()
					&& parent_set_ != parent_disjoint_set_.find_set(
							node_itr_->first));

			return *this;
		}

		int operator *() {
			return node_itr_->first;
		}

		bool operator==(const NodeIterator& other) {
			return (node_itr_ == other.node_itr_);
		}

		bool operator!=(const NodeIterator& other) {
			return (node_itr_ != other.node_itr_);
		}

	private:
		//PartIterator &part_iterator_;
		typename std::unordered_map<int, Element>::iterator node_itr_;
		int parent_set_;
		int num_nodes_seen_;
		int num_nodes_in_set_;
		DisjointSetForest &parent_disjoint_set_;
		//Node iterator needs to store either 1. reference to part itr
		// But will need to iterate through map without modifying part itr
		// Makes more sense to just store pointer to map
		// Find set
		// Then iterate based on that
	};

	//#################### PUBLIC METHODS ####################
public:
	/**
	 @brief  Adds a single element x (and its associated value) to the disjoint set forest.

	 @param[in]  x       The index of the element
	 @param[in]  value   The value to initially associate with the element
	 @pre
	 -   x must not already be in the disjoint set forest
	 */
	void add_element(int x, const T& value = T()) {
		m_elements.insert(std::make_pair(x, Element(value, x)));
		++m_setCount;
	}

	/**
	 @brief  Adds multiple elements (and their associated values) to the disjoint set forest.

	 @param[in]  elements    A map from the elements to add to their associated values
	 @pre
	 -   None of the elements to be added must already be in the disjoint set forest
	 */
	void add_elements(const std::map<int, T>& elements) {
		for (typename std::map<int, T>::const_iterator it = elements.begin(),
				iend = elements.end(); it != iend; ++it) {
			m_elements.insert(std::make_pair(it->first, Element(it->second,
					it->first)));
		}
		m_setCount += elements.size();
	}

	/**
	 @brief  Returns the number of elements in the disjoint set forest.

	 @return As described
	 */
	int element_count() const {
		return static_cast<int> (m_elements.size());
	}

	/**
	 @brief  Finds the index of the root element of the tree containing x in the disjoint set forest.

	 @param[in]  x   The element whose set to determine
	 @pre
	 -   x must be an element in the disjoint set forest
	 @throw Exception
	 -   If the precondition is violated
	 @return As described
	 */
	int find_set(int x) const {
		//std::cout << "last parent of 7 was " << get_parent(7) << std::endl;
		Element& element = get_element(x);
		int parent = element.m_parent;
		if (parent != x) {
			parent = find_set(parent);
		}
		//std::cout << "last parent of 7 was " << get_parent(7) << std::endl;

		return parent;
	}

	/**
	 @brief  Returns the current number of disjoint sets in the forest (i.e. the current number of trees).

	 @return As described
	 */
	int set_count() const {
		return m_setCount;
	}

	/**
	 @brief  Merges the disjoint sets containing elements x and y.

	 If both elements are already in the same disjoint set, this is a no-op.

	 @param[in]  x   The first element
	 @param[in]  y   The second element
	 hbyg    @pre
	 -   Both x and y must be elements in the disjoint set forest
	 @throw Exception
	 -   If the precondition is violated
	 */
	void union_sets(int x, int y) {
		int setX = find_set(x);
		int setY = find_set(y);
		if (setX != setY) {
			link(setX, setY);
			was_last_union_effective = true;
		} else {
			was_last_union_effective = false;
		}
	}

	/**
	 @brief Undoes the last union set, to avoid making copies
	 */
	void undo_last_union() {
		if (was_last_union_effective == true) {
			//int old = last_unioned_node_->m_parent;
			//int index;
			//for (typename std::map<int,Element>::iterator itr = m_elements.begin(); itr != m_elements.end(); ++itr) {
			//	if (last_unioned_node_ == &(itr->second)) {
			//		index = itr->first;
			//	}
			//}
			last_unioned_node_->m_parent = last_unioned_parent_;

			//std::cout << "changing node " << index << " from " << old << " to parent to " << last_unioned_parent_ << std::endl;
			//Element& elementX = get_element(last_);
			//int alpha;
			++m_setCount;
			was_last_union_effective = false; // Can't undo twice
		}
	}

	/**
	 @brief  Returns the value associated with element x.

	 @param[in]  x   The element whose value to return
	 @pre
	 -   x must be an element in the disjoint set forest
	 @throw Exception
	 -   If the precondition is violated
	 @return As described
	 */
	T& value_of(int x) {
		return get_element(x).m_value;
	}

	/**
	 @brief  Returns the value associated with element x.

	 @param[in]  x   The element whose value to return
	 @pre
	 -   x must be an element in the disjoint set forest
	 @throw Exception
	 -   If the precondition is violated
	 @return As described
	 */
	const T& value_of(int x) const {
		return get_element(x).m_value;
	}

	/**
	 @brief  Returns the size of a set.

	 @param[in]  part part iterator to set
	 */
	int set_size(PartIterator &partitr) {
		int size = 0;
		for (NodeIterator nitr = partitr.begin(); nitr != partitr.end(); ++nitr) {
			++size;
		}
		return size;
	}

	static std::set<int> intersection(PartIterator &partitr1,
			PartIterator &partitr2) {
		std::set<int> intersection_set;
		for (NodeIterator n1itr = partitr1.begin(); n1itr != partitr1.end(); ++n1itr) {
			for (NodeIterator n2itr = partitr2.begin(); n2itr != partitr2.end(); ++n2itr) {
				if (*n1itr == *n2itr) {
					intersection_set.insert(*n1itr);
				}
			}
		}

		return intersection_set;
	}

	std::map<int, int> set_to_node;

	/**
	 @brief  Assigns node to a set id (which needn't be another node)
	 */
	void add_node_to_set(int node_id, int set_id) {
		std::map<int, int>::iterator itr = set_to_node.find(set_id);
		if (itr == set_to_node.end()) {
			set_to_node[set_id] = node_id;
		} else {
			this->union_sets(itr->second, node_id);
		}
	}

	//#################### PRIVATE METHODS ####################
private:
	int get_parent(int x) const {
		Element &e = get_element(x);
		return (e.m_parent);
	}

	Element& get_element(int x) const {
		typename std::unordered_map<int, Element>::iterator it =
				m_elements.find(x);
		//TODO add exception handling
		if (it != m_elements.end())
			return it->second;
		else {
			std::cout << "Error nonexistent element: " << x << std::endl;
			return it->second;
		}
		//else throw Exception(OSSWrapper() << "No such element: " << x);
	}

	void link(int x, int y) {
		Element& elementX = get_element(x);
		Element& elementY = get_element(y);
		int& rankX = elementX.m_rank;
		int& rankY = elementY.m_rank;
		if (rankX > rankY) {
			last_unioned_node_ = &elementY;
			last_unioned_parent_ = elementY.m_parent;
			elementY.m_parent = x;
		} else {
			last_unioned_node_ = &elementX;
			last_unioned_parent_ = elementX.m_parent;
			elementX.m_parent = y;
			if (rankX == rankY)
				++rankY;
		}
		--m_setCount;
	}
};
}

#endif
