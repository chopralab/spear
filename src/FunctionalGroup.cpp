#include "spear/FunctionalGroup.hpp"
#include "spear/Molecule_impl.hpp"

#include <boost/graph/vf2_sub_graph_iso.hpp>

using namespace Spear;

static std::string str_toupper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), 
                   [](unsigned char c){ return std::toupper(c); }
                  );
    return s;
}

FunctionalGroup::FunctionalGroup(const std::string& smiles) {
    using boost::add_vertex;
    assert(!smiles.empty());

    std::queue<size_t> paren_backs;
    std::map <size_t, size_t> rings;
    size_t current = 0;
    size_t previous = 0;
    auto bond_order = Bond::Order::UNKNOWN;
    bool first_atom = true;
    bool in_prop_list = false;

    auto add_atom = 
    [this, &first_atom, &bond_order, &previous, &current](Element::Symbol atom) {
        add_vertex(atom, graph_);

        if (!first_atom) {
            boost::add_edge(previous, ++current,
            EdgeProperty(bond_order), graph_);
        }

        first_atom = false;
        previous = current;
        bond_order = Bond::Order::UNKNOWN;
        properties_.push_back(std::list<AtomPropertyCompare>());
    };

    auto read_number =
    [&smiles](size_t& start) {
        if (!std::isdigit(smiles[start + 1])) {
            return static_cast<size_t>(1);
        }
        ++start;
        size_t number = 0;
        while (start < smiles.size() && std::isdigit(smiles[start])) {
            number *= 10;
            number += static_cast<size_t>(smiles[start] - '0');
            ++start;
        }
        --start;
        return number;
    };

    for (size_t i = 0; i < smiles.size(); ++i) {
        if (in_prop_list) {
            size_t bonds = 0;
            std::cout << smiles[i] << std::endl;
            switch (smiles[i]) {
            case 'D':
                bonds = read_number(i);

                // Add a lambda function to compare the bond counts
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        return a1.degree() == bonds;
                    }
                );
                continue;
            case 'X':
                bonds = read_number(i);
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        return a1.expected_bonds() == bonds;
                    }
                );
                continue;
            case 'h':
                bonds = read_number(i);
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        return a1.implicit_hydrogens() == bonds;
                    }
                );
                continue;
            case 'H':
                bonds = read_number(i);
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        return a1.explicit_hydrogens() == bonds;
                    }
                );
                continue;
            case 'R':
                bonds = read_number(i);
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        auto sssrs = a1.sssrs();
                        return std::distance(sssrs.first, sssrs.second) == bonds;
                    }
                );
                continue;
            case 'r':
                bonds = read_number(i);
                properties_.back().emplace_back(
                    [bonds](const AtomVertex& a1) {
                        auto sssrs = a1.sssrs();
                        for (auto ring = sssrs.first; ring != sssrs.second; ++ring) {
                            if (ring->second.size() == bonds) return true;
                        }
                        return false;
                    }
                );
                continue;
            case 'a':
                properties_.back().emplace_back(
                    [](const AtomVertex& a1) {
                        return a1.is_aromatic();
                    }
                );
                continue;
            case 'A':
                properties_.back().emplace_back(
                    [](const AtomVertex& a1) {
                        return !a1.is_aromatic();
                    }
                );
                continue;
            default:
                break;
            }
        }

        if (smiles[i] == 'a' || smiles[i] == 'A') {
            bool should_be_aromatic = std::islower(smiles[i]);
            if (!in_prop_list) {
                add_atom(Element::Symbol(0));
            }

            properties_.back().emplace_back(
                [should_be_aromatic](const AtomVertex& a1) {
                    return a1.is_aromatic() == should_be_aromatic;
                }
            );

            continue;
        }

        if (std::islower(smiles[i])) {
            auto is_found = Element::SymbolForName.find(
                str_toupper(smiles.substr(i, 1))
            );
            if (is_found != Element::SymbolForName.end()) {
                add_atom(is_found->second);
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
            auto is_found = Element::SymbolForName.find(element_name);

            if (is_found != Element::SymbolForName.end()) {
                add_atom(is_found->second);
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
        case '*': add_atom(Element::Symbol(0)); break;
        case '@': break;
        case '/': break;
        case '\\': break;
        case '%': break;
        case '~': bond_order = Bond::Order::UNKNOWN; break;
        case '-': bond_order = Bond::Order::SINGLE; break; // We don't support charges
        case '=': bond_order = Bond::Order::DOUBLE; break;
        case '#': bond_order = Bond::Order::TRIPLE; break;
        case ':': bond_order = Bond::Order::AROMATIC; break;
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

std::list<std::vector<size_t>> Spear::find_functional_groups(const Molecule& mol,
                                                             const FunctionalGroup& fg) {

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
