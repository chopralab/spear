// License:

#ifndef SPEAR_FUNCTIONALGROUP_HPP
#define SPEAR_FUNCTIONALGROUP_HPP

#include <queue>
#include <map>

#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "spear/exports.hpp"
#include "spear/Molecule.hpp"
#include "chemfiles/periodic_table.hpp"
#include "chemfiles/Connectivity.hpp"

namespace Spear {

using AtomPropertyCompare = std::function<bool(const AtomVertex&)>;

class FunctionalGroup {
    Graph graph_;
    std::list<std::list<AtomPropertyCompare>> properties_;
public:
    FunctionalGroup(const std::string& smiles);

    const Graph& graph() const {
        return graph_;
    }

    const std::list<AtomPropertyCompare>& properties(size_t i) const {
        auto start = properties_.cbegin();
        std::advance(start, i);
        return *start;
    }
};

std::list<std::vector<size_t>> find_functional_groups(const Molecule& mol,
                                                      const FunctionalGroup& fg);

}

#endif
