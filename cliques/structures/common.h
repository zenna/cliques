#pragma once

/*
 TODO: this common should not depend on partition,
 use templated typedefs instead http://stackoverflow.com/questions/26151/template-typedefs-whats-your-work-around
 or get rid of typedefs alltogether
 */
namespace cliques {

template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

/**
 @brief  Hashing functor for general partition

 Functor creates a hash which should work for all partitions
 It is probably wise to use template specialisation to write new hashers
 for different partition implementation
 */
struct partition_hash {
    template<typename P>
    size_t operator()(P &partition) const {
        std::hash<int> ihash;
        std::size_t seed = 0;
        int num_elements = partition.element_count();
        for (int i = 0; i < num_elements; ++i) {
            hash_combine(seed, partition.find_set(i));
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