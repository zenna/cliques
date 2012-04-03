#pragma once

#if defined USE_BOOST
#include <boost/functional/hash.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#else
#include <unordered_map>
#include <unordered_set>
#endif

/*
 TODO: this common should not depend on partition,
 use templated typedefs instead http://stackoverflow.com/questions/26151/template-typedefs-whats-your-work-around
 or get rid of typedefs alltogether
 */
namespace cliques {

template<typename T>
struct xyz_colour {
    xyz_colour(T x, T y, T z) :
        x(x), y(y), z(y) {
    }
    xyz_colour() :
        x(0), y(0), z(0), a(0) {
    }
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
struct partition_hash {
    template<typename P>
    size_t operator()(P &partition) const {
        boost::hash<int> ihash;
        std::size_t seed = 0;
        int num_elements = partition.element_count();
        for (int i = 0; i < num_elements; ++i) {
            boost::hash_combine(seed, partition.find_set(i));
        }
        return ihash(seed);
    }
};

/**
 @brief  Equality functor for general partition

 */
struct partition_equal {
    template<typename P>
    bool operator()(P const& x, P const& y) const {
        return (x == y);
    }
};

}
