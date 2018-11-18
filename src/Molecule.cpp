#include "spear/Molecule.hpp"

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/hawick_circuits.hpp>
#include <boost/graph/undirected_graph.hpp>

using namespace Spear;

// The cycle_printer is a visitor that will print the path that comprises
// the cycle. Note that the back() vertex of the path is not the same as
// the front(). It is implicit in the listing of vertices that the back()
// vertex is connected to the front().
struct cycle_saver
{

    cycle_saver(std::set<std::set<size_t>>& counter):
        found_rings(counter){
    }

    template <typename Path, typename Graph>
    void cycle(const Path& p, const Graph& g) {
        // Get the property map containing the vertex indices that are saved
        auto indices = boost::get(boost::vertex_index, g);

        // We deal mostly with undirected graphs, so a bond is technically a
        // cycle. Let's avoid those!
        if (p.size() <= 2) {
            return;
        }

        // We've found a ring, so add it to found rings
        std::set<size_t> new_ring;
        for(const auto& i : p) {
            new_ring.insert(boost::get(indices, i));
        }

        found_rings.insert(std::move(new_ring));
    }

    std::set<std::set<size_t>>& found_rings;
};

Molecule::Molecule(const chemfiles::Frame& frame) :
    frame_(frame.clone()), graph_() {

    auto &topo = frame_.topology();

    std::map<size_t, VertexDescriptor> vertices;
    for (auto atom : frame) {
        boost::add_vertex(VertexDescriptor{*(atom.atomic_number())}, graph_);
    }

    // Create the graph representation
    auto bonds = topo.bonds();
    for ( size_t i = 0; i < bonds.size(); ++i) {
        boost::add_edge(bonds[i][0], bonds[i][1],
            EdgeProperty(topo.bond_orders()[i]), graph_);
    }
}

const std::set<std::set<size_t>> Molecule::rings() const {
    std::set<std::set<size_t>> ret_rings;

    cycle_saver vis(ret_rings);
    boost::hawick_circuits(graph_, vis);

    return ret_rings;
}

std::vector<Spear::EdgeDescriptor> Molecule::get_bonds_in(const std::set<size_t>& atoms) const {
    EdgeIterator begin, end;
    std::tie(begin, end) = boost::edges(graph_);

    std::vector<EdgeDescriptor> ret;

    std::copy_if(begin, end, std::back_inserter(ret),
        [&atoms, this](EdgeDescriptor e){
            auto target = boost::target(e, graph_);
            auto source = boost::source(e, graph_);
            if (atoms.count(target) == 0) {
                return false;
            }
            if (atoms.count(source) == 0) {
                return false;
            }
            return true;
    });

    return ret;
}
