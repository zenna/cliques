#ifndef CLIQUES_BITSETPARTITION_H
#define CLIQUES_BITSETPARTITION_H

#include <cmath>
#include <map>
#include <set>
#include <iostream>

#include <boost/dynamic_bitset.hpp>

namespace cliques {

/**
@brief  A space efficient partition using a boost dynamic bitset

This partition class is designed for space efficiency, internally
it's structure is a simple binary encoded list of numbers in the form:
e.g. 01230 means that node 0 and node 5 are in the same group,
whereas all other nodes are in their own singleton groups
*/
class BitsetPartition
    {
	//Idea partition would
	   // minimal storage requirements less than 8 bytes
	   // Fast iteration through parts and nodes
	   // Fast
        //#################### PRIVATE VARIABLES ####################
    private:
        boost::dynamic_bitset<> partition;
    	int bits_per_block;
    	int num_nodes;


        //#################### CONSTRUCTORS ####################
    public:
        /**
        @brief  Constructs an partition with each node in its own group.
        */
        BitsetPartition(int num_nodes) : num_nodes(num_nodes) {
        	bits_per_block = ceil(float(std::log(num_nodes) / std::log(2)));
        	partition = boost::dynamic_bitset(num_bits*num_nodes);
        }

        //#################### ITERATORS ####################
        class NodeIterator;

        class PartIterator {
        public:
            PartIterator(int current_set, BitsetPartition &parent_bitset_partition)
                : current_set(current_set), parent_bitset_partition(parent_bitset_partition) {}

            PartIterator& operator=(const PartIterator& other)
            {
            	current_set = other.current_set;
                visited_sets = other.visited_sets;
                return(*this);
            }

            bool operator==(const PartIterator& other)
            {
                return(current_set == other.current_set);
            }

            bool operator!=(const PartIterator& other)
            {
                return(part_itr_ != other.part_itr_);
            }

            PartIterator& operator++() {
            	// Increase unless I've seen before

            	int parent_set = parent_disjoint_set_.find_set(part_itr_->first);
                visited_parents_.insert(parent_set);

                if (visited_parents_.size() == parent_disjoint_set_.m_setCount) {
                    part_itr_ = parent_disjoint_set_.m_elements.end();
                }

                else {
                    do
                    {
                        ++part_itr_;
                        parent_set = parent_disjoint_set_.find_set(part_itr_->first);
                    }
                    while(visited_parents_.find(parent_set) != visited_parents_.end());
                    //visited_parents_.insert(parent_set);
                }

                return *this;
            }

            int operator*()
            {
                return(parent_disjoint_set_.find_set(part_itr_->first));
            }

            NodeIterator begin() {
                typename std::map<int,Element>::iterator tmp_node_itr_ = part_itr_;
                return NodeIterator(tmp_node_itr_ , this->parent_disjoint_set_);
            }

            NodeIterator end() {
                typename std::map<int,Element>::iterator tmp_node_itr_ = parent_disjoint_set_.m_elements.end();
                return NodeIterator(tmp_node_itr_ , this->parent_disjoint_set_);
            }

        private:
            //typename std::map<int,Element>::iterator part_itr_;
            int current_set;
            //typename std::map<int,T>::iterator part_itr_;
            //std::set<int> visited_parents_;
            BitsetPartition &parent_bitset_partition;
            //Node* node_;
        };

        PartIterator begin() {
        	node_position = 0;
            return PartIterator(node_position, *this);
        }

        PartIterator end() {
        	node_position = this->partition.size();
        	return PartIterator(node_position, *this);
        }

        class NodeIterator {
        public:
            NodeIterator(typename std::map<int,Element>::iterator node_itr, DisjointSetForest &parent_disjoint_set)
                : node_itr_(node_itr), parent_disjoint_set_(parent_disjoint_set) {
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
                do
                {
                    ++node_itr_;
                }
                while(node_itr_ != parent_disjoint_set_.m_elements.end() &&
                		parent_set_ != parent_disjoint_set_.find_set(node_itr_->first));

                return *this;
            }

            int operator *() {
                return node_itr_->first;
            }

            bool operator==(const NodeIterator& other)
            {
                return(node_itr_ == other.node_itr_);
            }

            bool operator!=(const NodeIterator& other)
            {
                return(node_itr_ != other.node_itr_);
            }

        private:
            //PartIterator &part_iterator_;
            typename std::map<int,Element>::iterator node_itr_;
            int current_set;
            int num_nodes_seen_;
            int num_nodes_in_set_;
            DisjointSetForest &parent_disjoint_set_;
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
        void add_element(int x)
        {
            m_elements.insert(std::make_pair(x, Element(value, x)));
            ++m_setCount;
        }

        /**
        @brief  Returns the number of elements in the disjoint set forest.

        @return As described
        */
        int element_count() const
        {
            return num_nodes;
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
        int find_set(int x) const
        {
        	boost::dynamic_bitset<> node_set(bits_per_node);
        	for (int i=0;i<bits_per_node;++i) {
        		node_set[i] = partition[bits_per_node * node_id + i];
        	}
        	return node_set.to_ulong();
        }

        /**
        @brief  Returns the current number of disjoint sets in the forest (i.e. the current number of trees).

        @return As described
        */
        int set_count() const
        {
            return 5;
        }

        /**
        @brief  Merges the disjoint sets containing elements x and y.

        If both elements are already in the same disjoint set, this is a no-op.

        @param[in]  x   The first element
        @param[in]  y   The second element
        */
        void union_sets(int node_x, int node_y)
        {
        	// pass through list, change one to other
        	boost::dynamic_bitset<> set_of_x = find_bitset(node_x);
        	boost::dynamic_bitset<> set_of_y = find_bitset(node_y);

        	for (int i=0;i<num_nodes; ++i) {
        		if (find_bitset(i) == set_of_x) {
        			add_node_to_set(i, set_of_y);
        		}
        	}
        }
        /**
        @brief  Returns the size of a set.

        @param[in]  part part iterator to set
        */
        int set_size(PartIterator &partitr)
        {
			int size = 0;
			for (NodeIterator nitr = partitr.begin(); nitr != partitr.end(); ++nitr) {
				++size;
			}
			return size;
        }

        static std::set<int> intersection(PartIterator &partitr1, PartIterator &partitr2) {
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

        /**
        @brief  Assigns node to a set id (which needn't be another node)
        */
        void add_node_to_set(int node_id, int set_id) {
        	const boost::dynamic_bitset<> set_id_bitset(bits_per_block, set_id);
        	for (int i=0;i<bits_per_block;++i) {
        		partition[bits_per_block * node_id + i] = set_id_bitset[i];
        	}
        }

        //#################### PRIVATE METHODS ####################
    private:
        boost::dynamic_bitset<> find_bitset(int node_id) {
        	boost::dynamic_bitset<> node_set(bits_per_node);
        	for (int i=0;i<bits_per_node;++i) {
        		node_set[i] = partition[bits_per_node * node_id + i];
        	}
        	return node_set;
        }

        void add_node_to_bitset(int node_id, boost::dynamic_bitset &set) {
        	for (int i=0;i<bits_per_block;++i) {
        		partition[bits_per_block * node_id + i] = set[i];
        	}
        }

    };
}

#endif
