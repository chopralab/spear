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

class FunctionalGroup {
    Graph graph_;

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

        for (size_t i = 0; i < smiles.size(); ++i) {
            if (std::isupper(smiles[i]) || smiles[i] == '*') {
                if (smiles[i] != '*') {
                    size_t element_length =
                        i + 1 < smiles.size() && std::islower(smiles[i+1])? 2 : 1;

                    auto element_name = smiles.substr(i, element_length);

                    auto temp = chemfiles::Atom("", element_name);
                    if (temp.atomic_number()) {
                        add_vertex(*(temp.atomic_number()), graph_);
                    } else {
                        throw std::invalid_argument("Element not found: " + element_name);
                    }
                    i += element_length - 1;
                } else {
                    add_vertex(0, graph_);
                }

                if (!first_atom) {
                    boost::add_edge(previous, ++current,
                        EdgeProperty(bond_order), graph_);
                }

                first_atom = false;
                previous = current;
                bond_order = Bond::UNKNOWN;
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
            case '.': break;
            case '@': break;
            case '/': break;
            case '\\': break;
            case '%': break;
            case '-': bond_order = Bond::SINGLE; break;
            case '=': bond_order = Bond::DOUBLE; break;
            case '#': bond_order = Bond::TRIPLE; break;
            case ':': bond_order = Bond::AROMATIC; break;
            case '[': break;
            case ']': break;
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
};

struct FunctionalGroupFinder {
    FunctionalGroupFinder(const Molecule& mol, const FunctionalGroup& fg,
        std::list<std::vector<size_t>>& found_groups):
    mol_(mol), fg_(fg), found_groups_(found_groups) {
    }

    template <typename Map1To2, typename Map2To1>
    bool operator()(Map1To2 f, Map2To1) {

        std::vector<size_t> group;
        group.reserve(boost::num_vertices(fg_.graph()));
        auto verticies_iter = boost::vertices(fg_.graph());

        for (auto v = verticies_iter.first; v != verticies_iter.second; ++v) {
            group.push_back(boost::get(boost::vertex_index_t(), mol_.graph(), boost::get(f, *v)));
        }

        found_groups_.emplace_back(group);
        return true;
    }

private:
    const Molecule& mol_;
    const FunctionalGroup& fg_;
    std::list<std::vector<size_t>>& found_groups_;
};

template <typename PropertyFirst, typename PropertySecond, typename IgnoreType>
struct CompareProperty {
  
    CompareProperty(const PropertyFirst property_map1,
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

template <typename PropertyFirst, typename PropertySecond, typename IgnoreType>
CompareProperty<PropertyFirst, PropertySecond, IgnoreType>
compare_property(const PropertyFirst property1,
                 const PropertySecond property2,
                 const IgnoreType ignore) {
    return CompareProperty<PropertyFirst, PropertySecond, IgnoreType>
        (property1, property2, ignore);
}

std::list<std::vector<size_t>> find_functional_groups(const Molecule& mol, const FunctionalGroup& fg) {
    std::list<std::vector<size_t>> found_groups;
    FunctionalGroupFinder callback(mol, fg, found_groups);

    const auto& graph1 = fg.graph();
    const auto& graph2 = mol.graph();

    auto vertex_comp =
        compare_property(boost::get(boost::vertex_name, graph1),
                         boost::get(boost::vertex_name, graph2),
                         uint64_t(0));

    auto edge_comp =
        compare_property(boost::get(boost::edge_name, graph1),
                         boost::get(boost::edge_name, graph2),
                         uint64_t(0));

    boost::vf2_subgraph_iso(graph1, graph2, callback,
        boost::vertex_order_by_mult(graph1),
        boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));

    return found_groups;
}

}

#endif
