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
        size_t bond_order = Bond::SINGLE;

        auto temp = chemfiles::Atom("", smiles.substr(0,1));
        add_vertex(*(temp.atomic_number()), graph_);

        for (size_t i = 1; i < smiles.size(); ++i) {
            if (std::isupper(smiles[i])) {
                ++current;
                temp = chemfiles::Atom("", smiles.substr(i,1));
                add_vertex(*(temp.atomic_number()), graph_);

                boost::add_edge(previous, current,
                    EdgeProperty(bond_order), graph_);

                previous = current;
                bond_order = Bond::SINGLE;
                continue;
            }

            if (std::isdigit(smiles[i])) {
                size_t ring_id = static_cast<size_t>(smiles[i]);
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

std::list<std::vector<size_t>> find_functional_groups(const Molecule& mol, const FunctionalGroup& fg) {
    std::list<std::vector<size_t>> found_groups;
    FunctionalGroupFinder callback(mol, fg, found_groups);

    const auto& graph1 = fg.graph();
    const auto& graph2 = mol.graph();

    auto vertex_comp =
    boost::make_property_map_equivalent(boost::get(boost::vertex_name, graph1),
                                        boost::get(boost::vertex_name, graph2));

    auto edge_comp =
    boost::make_property_map_equivalent(boost::get(boost::edge_name, graph1),
                                        boost::get(boost::edge_name, graph2));

    boost::vf2_subgraph_iso(graph1, graph2, callback,
        boost::vertex_order_by_mult(graph1),
        boost::edges_equivalent(edge_comp).vertices_equivalent(vertex_comp));

    return found_groups;
}

}

#endif
