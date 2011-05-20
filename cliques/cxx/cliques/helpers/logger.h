#pragma once

namespace cliques {

struct NoLogging {
    template <typename P>
    void log(const P &) {}
};

template <typename P>
struct Logging {
    std::vector<P> vec;

    //template <typename P>
    void log(const P &p) {
        vec.push_back(p);
    }
};

}
