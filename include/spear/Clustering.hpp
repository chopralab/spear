// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_CLUSTERING_HPP
#define SPEAR_CLUSTERING_HPP

#include <vector>
#include <algorithm>

namespace Spear {

template<typename Container, typename DistanceFunction>
void greedy_remove_non_representative(Container& conformations,
                                      DistanceFunction func, double radius) {
    for (size_t i = 1; i < conformations.size(); ++i) {
        // i starts at 1, so next is the first configuration that can be removed!
        auto next = conformations.begin() + i;
        conformations.erase(std::remove_if(next, conformations.end(),
            [&func, radius, next, &conformations](const typename Container::reference current) {
                return radius > func(*(next - 1), current);
        }), conformations.end());
    }
}

}

#endif
