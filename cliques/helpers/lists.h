/* Copyright (c) Z Tavares, M Schaub - zennatavares@gmail.com, 2010-2012
Helper functions for working with lists, vectors and other array like containers */
#pragma once
namespace cliques {

/**
 @brief  Summation for iterable containers of numerical type
 @tparam cN                         Numerical type (that can be summed)
 @param[in]  container              container containing the values, e.g. a vector of doubles
 @param[out] total                  The sum (i.e. total)
 */
template<typename C>
typename C::value_type sum(C container ) {
    typedef typename C::value_type value_type;
    // Value_type of container
    value_type total;
    for (value_type value : container) {
        total += value;
    }
    return total;
}
}
