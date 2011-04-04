#ifndef CLIQUES_COMMON_H
#define CLIQUES_COMMON_H

#include <bitset>
#include <set>
#include <structures/partition.h>

#if defined USE_BOOST
    #include <boost/functional/hash.hpp>
    #include <boost/unordered_map.hpp>
    #include <boost/unordered_set.hpp>
#else
    #include <unordered_map>
    #include <unordered_set>
#endif

#define N 128
#define MAX_LINE_LENGTH 1000
#define SET_MAX_SIZE 63

#define NUM_NODES 16
#define BITS_PER_NODE 4

/*
TODO: this common should not depend on partition,
use templated typedefs instead http://stackoverflow.com/questions/26151/template-typedefs-whats-your-work-around
or get rid of typedefs alltogether
*/
namespace cliques {

#if defined USE_BOOST
	typedef boost::hash<std::string> shash;
#else
	typedef std::hash<std::string> shash;
#endif
struct bitset_hash
{
	size_t operator()(std::bitset<N> bitset) const
	{
		return shash()(bitset.to_string());
	}
};

#if defined USE_BOOST
	typedef boost::unordered_map<std::bitset<N>, cliques::Partition, bitset_hash> umap;
#else
	typedef std::unordered_map<std::bitset<N>, cliques::Partition, bitset_hash> umap;
#endif

struct iterator_hash
{
	size_t operator()(umap::const_iterator it) const
	{
		return shash()(it->first.to_string());
	}
};

#if defined USE_BOOST
	typedef boost::unordered_set<umap::const_iterator, iterator_hash> uset;
#else
	typedef std::unordered_set<umap::const_iterator, iterator_hash> uset;
#endif

struct set_hash
{
	size_t operator()(std::set<int> set_key) const
	{
		int alpha = 0;
		for (std::set<int>::iterator itr = set_key.begin(); itr != set_key.end(); ++itr)
		  alpha += *itr;
		return alpha; //shash()(bitset.to_string());
	}
};

#if defined USE_BOOST
	typedef boost::unordered_map<std::set<int>,float,set_hash> setumap;
#else
	typedef std::unordered_map<std::set<int>,float,set_hash> setumap;
#endif

template <typename T>
struct xyz_colour {
	xyz_colour(T x, T y, T z) : x(x), y(y), z(y) {}
	xyz_colour() : x(0), y(0), z(0), a(0) {}
	T x;
	T y;
	T z;
	T a;
};

/**
@brief  Hashing functor for general partition

Functor creates a hash which should work for all partitions
It is probably wise to use template specialisation to write new hashers
for different partition implementation
*/
struct partition_hash
{
	template <typename P>
	size_t operator()(P &partition) const
	{
		boost::hash<int> ihash;
		return ihash(partition.set_count());
	}
};

/**
@brief  Equality functor for general partition

*/
struct partition_equal
{
    template <typename P>
    bool operator()(P const& x,
        P const& y) const
    {
        return (x.set_to_node == y.set_to_node);
    }
};


} //namespace cliques

#endif //CLIQUES_STRUCTURES_H