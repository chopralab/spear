// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_RINGS_HPP
#define SPEAR_RINGS_HPP

#include <set>
#include <unordered_map>

namespace Spear {

struct ring_compare_functor {
    bool operator()(const std::set<std::size_t>& a,
                    const std::set<std::size_t>& b) const {
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

using AtomRingMap = std::unordered_multimap<size_t, const std::set<size_t>&>;
using AtomRingMapIterator = AtomRingMap::const_iterator;
using AtomRingMapIteratorPair = std::pair<AtomRingMapIterator,AtomRingMapIterator>;

}

#endif
