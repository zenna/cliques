#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace cliques {
std::vector<char> string_to_char_p(std::string str) {
    std::vector<char> writable(str.size() + 1);
    std::copy(str.begin(), str.end(), writable.begin());
    return writable;
}

template<class T>
inline std::string to_string(const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}

template<class T>
inline std::string to_char_pointer(const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}

}
