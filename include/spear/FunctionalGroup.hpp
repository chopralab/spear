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
    FunctionalGroup(const std::string& smiles) {
        using chemfiles::find_in_periodic_table;
        using chemfiles::Bond;
        using boost::add_vertex;
        assert(!smiles.empty());

        std::queue<size_t> paren_backs;
        std::map <size_t, size_t> rings;
        size_t current = 0;
        size_t previous = 0;
        size_t bond_order = Bond::UNKNOWN;
        bool first_atom = true;
        bool in_prop_list = false;

        auto add_atom = 
        [this, &first_atom, &bond_order, &previous, &current](uint64_t atom) {
            add_vertex(atom, graph_);

            if (!first_atom) {
                boost::add_edge(previous, ++current,
                EdgeProperty(bond_order), graph_);
            }

            first_atom = false;
            previous = current;
            bond_order = Bond::UNKNOWN;
            properties_.push_back(std::list<AtomPropertyCompare>());

        };

        for (size_t i = 0; i < smiles.size(); ++i) {
            if (in_prop_list && smiles[i] == 'D' && i != smiles.size() - 1) {
                auto bonds = static_cast<size_t>(smiles[i + 1] - '0');

                // Add a lambda function to compare the bond counts
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        return a1.neighbor_count() == bonds;
                    }
                );
                i++;
                continue;
            }

            if (smiles[i] == 'a' || smiles[i] == 'A') {
                bool should_be_aromatic = std::islower(smiles[i]);
                if (!in_prop_list) {
                    add_atom(0);
                }

                properties_.back().emplace_back(
                    [should_be_aromatic](const AtomVertex& a1) {
                        return a1.is_aromatic() == should_be_aromatic;
                    }
                );

                continue;
            }

            if (std::islower(smiles[i])) {
                auto temp = chemfiles::Atom("", smiles.substr(i, 1));
                if (temp.atomic_number()) {
                    add_atom(*(temp.atomic_number()));
                } else {
                    throw std::invalid_argument("Element not found: " + smiles[i]);
                }
                properties_.back().emplace_back(
                    [](const AtomVertex& a1) {
                        return a1.is_aromatic();
                    }
                );
                continue;
            }

            // Unlike a smiles string, C does not imply aliphatic!
            if (std::isupper(smiles[i])) {
                size_t element_length =
                    i + 1 < smiles.size() && std::islower(smiles[i+1])? 2 : 1;

                auto element_name = smiles.substr(i, element_length);

                auto temp = chemfiles::Atom("", element_name);
                if (temp.atomic_number()) {
                    add_atom(*(temp.atomic_number()));
                } else {
                    throw std::invalid_argument("Element not found: " + element_name);
                }
                i += element_length - 1;
                continue;
            }

            if (std::isdigit(smiles[i])) {
                size_t ring_id = static_cast<size_t>(smiles[i] - '0');
                auto ring_lookup = rings.find(ring_id);

                if (ring_lookup == rings.end()) {
                    rings.insert({ring_id, previous});
                    continue;
                }

                boost::add_edge(previous, ring_lookup->second,
                    EdgeProperty(bond_order), graph_);
                rings.erase(ring_lookup);
                continue;
            }

            switch (smiles[i]) {
            case '.': first_atom = false; break;
            case '*': add_atom(0); break;
            case '@': break;
            case '/': break;
            case '\\': break;
            case '%': break;
            case '~': bond_order = Bond::UNKNOWN; break;
            case '-': bond_order = Bond::SINGLE; break; // We don't support charges
            case '=': bond_order = Bond::DOUBLE; break;
            case '#': bond_order = Bond::TRIPLE; break;
            case ':': bond_order = Bond::AROMATIC; break;
            case '[': in_prop_list = true; break;
            case ']': in_prop_list = false; break;
            case '+': break;
            case '(':
                paren_backs.push(previous);
                break;
            case ')':
                previous = paren_backs.back();
                paren_backs.pop();
                break;
            default: break;
            }
        }
    }

    const Graph& graph() const {
        return graph_;
    }

    const std::list<AtomPropertyCompare>& properties(size_t i) const {
        auto start = properties_.cbegin();
        std::advance(start, i);
        return *start;
    }
};

/// This Functor is run when the graph matching algorithm finds a potential
/// FunctionalGroup. It is responsible for populating the found_groups list
/// and performing final checks for properties in the functional group.
struct FunctionalGroupFinder {
    FunctionalGroupFinder(const Molecule& mol, const FunctionalGroup& fg,
        std::list<std::vector<size_t>>& found_groups):
    mol_(mol), fg_(fg), found_groups_(found_groups) {
    }

    // Map1 should be the functional group
    // Map2 should be the molecule
    // Therefore, fg_to_mol coverts the functional group mapping to the molecule
    template <typename Map1To2, typename Map2To1>
    bool operator()(Map1To2 fg_to_mol, Map2To1) {

        std::vector<size_t> group;

        // Prepare to loop through the functional group
        group.reserve(boost::num_vertices(fg_.graph()));
        auto verticies_iter = boost::vertices(fg_.graph());

        // Add the **molecule** index to the group, *not* the functional groups
        // v is an interator through the functional group
        for (auto v = verticies_iter.first; v != verticies_iter.second; ++v) {

            const auto& mol_id = boost::get(fg_to_mol, *v); // v is the functional group index
            const auto& mol_index = boost::get(boost::vertex_index_t(), mol_.graph(), mol_id);

            // Final property check
            auto fg_index = boost::get(boost::vertex_index_t(), fg_.graph(), *v);
            const auto& props = fg_.properties(fg_index);
            for (const auto& prop : props) {
                if (!prop(mol_[mol_index])) {
                    // Don't add the group, get out!
                    return true;
                }
            }

            // All tests for this functional group atom -> molecular atom passed
            group.push_back(mol_index);
        }

        // All tests for all atoms have passed
        found_groups_.emplace_back(group);
        return true;
    }

private:
    const Molecule& mol_;
    const FunctionalGroup& fg_;
    std::list<std::vector<size_t>>& found_groups_;
};


/// This Functor compares the properties of an edge (bond) or vertex(atom).
/// The template arguments should be automatically deduced when a comparison
/// structure is initialized using the **compare_topology** function.
template <typename PropertyFirst, typename PropertySecond, typename IgnoreType>
struct CompareTopology {
  
    CompareTopology(const PropertyFirst property_map1,
                    const PropertySecond property_map2,
                    const IgnoreType ignore) :
    property1_(property_map1), property2_(property_map2), ignore_(ignore) {}

    template <typename ItemFirst,
              typename ItemSecond>
    bool operator()(const ItemFirst item1, const ItemSecond item2) const {
        if (boost::get(property1_, item1) == ignore_) return true;
        if (boost::get(property2_, item2) == ignore_) return true;
        return boost::get(property1_, item1) == boost::get(property2_, item2);
    }
  
private:
    const PropertyFirst property1_;
    const PropertySecond property2_;
    const IgnoreType ignore_;
};

/// This is a convience function for creating CompareTopology Functors.
template <typename PropertyFirst, typename PropertySecond, typename IgnoreType>
CompareTopology<PropertyFirst, PropertySecond, IgnoreType>
compare_topology(const PropertyFirst property1,
                 const PropertySecond property2,
                 const IgnoreType ignore) {
    return CompareTopology<PropertyFirst, PropertySecond, IgnoreType>
        (property1, property2, ignore);
}

std::list<std::vector<size_t>> find_functional_groups(const Molecule& mol, const FunctionalGroup& fg) {

    // Initial objects to hold the found functional groups
    std::list<std::vector<size_t>> found_groups;
    FunctionalGroupFinder callback(mol, fg, found_groups);

    const auto& fg_graph = fg.graph();
    const auto& mol_graph = mol.graph();

    // Make comparison functions
    auto vertex_comp =
        compare_topology(boost::get(boost::vertex_name, fg_graph),
                         boost::get(boost::vertex_name, mol_graph),
                         uint64_t(0));

    auto edge_comp =
        compare_topology(boost::get(boost::edge_name, fg_graph),
                         boost::get(boost::edge_name, mol_graph),
                         uint64_t(0));

    // Find subgraphs
    boost::vf2_subgraph_iso(fg_graph, mol_graph, callback,
        boost::vertex_order_by_mult(fg_graph),
        boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));

    return found_groups;
}

}

#endif
