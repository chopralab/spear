#ifndef SPEAR_RINGS_HPP
#define SPEAR_RINGS_HPP

#include <set>

namespace Spear {

struct ring_compare_functor {
    bool operator()(const std::set<std::size_t>& a,
                    const std::set<std::size_t>& b) {
        if (a.size() < b.size()) {
            return true;
        }

        if (b.size() < a.size()) {
            return false;
        }

        return a < b;
    }
};

using RingSet = std::set<std::set<std::size_t>, ring_compare_functor>;

}

#endif
