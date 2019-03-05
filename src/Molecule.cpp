#include "spear/Molecule.hpp"

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/hawick_circuits.hpp>
#include <boost/graph/undirected_graph.hpp>

#include "spear/AtomType.hpp"

using namespace Spear;

size_t AtomVertex::expected_bonds() const {
    auto atomic_num = atomic_number();

    switch(atomic_num) {
        case 1:  // Hydrogen
        case 9:  // Fluorine
        case 17: // Chlorine
        case 35: // Bromine
        case 53: // Iodine
        return 1;
        break;

        case 5:  // Boron
        case 13: // Aluminium
        return 3;
        break;

        default: break;
    }

    auto hybridization_state = br_->get_default_atomtype()->hybridization(index_);

    if(atomic_num == 6) { // Carbon
        switch (hybridization_state) {
            case Hybridization::SP: return 2; break;
            case Hybridization::SP2: return 3; break;
            case Hybridization::SP3: return 4; break;
            default: return 0; break;
        }
    }

    if(atomic_num == 7) { // Nitrogen
        switch (hybridization_state) {
            case Hybridization::SP: return 1; break;
            case Hybridization::SP2: return 2; break;
            case Hybridization::SP3: return 3; break;
            default: return 0; break;
        }
    }

    if(atomic_num == 8) { // Oxygen
        switch (hybridization_state) {
            case Hybridization::SP: return 1; break;
            case Hybridization::SP2: return 1; break;
            case Hybridization::SP3: return 2; break;
            default: return 0; break;
        }
    }

    if(atomic_num == 15 || atomic_num == 16) { // Phosphorus and Sulfur
        if (is_aromatic()) return 2;
        return neighbor_count(); // Hard to tell, but things shouldn't be protonated
    }

    // Give up....
    return static_cast<size_t>(hybridization_state);
}

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

void Molecule::init_() {
    auto &topo = frame_.topology();

    std::map<size_t, VertexDescriptor> vertices;
    for (auto atom : frame_) {
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

static double cos_sim(const chemfiles::Vector3D& u, const chemfiles::Vector3D& v) {
    auto arc = chemfiles::dot(u,v) / (u.norm() * v.norm());
    return arc >= 1? 0 : std::acos(arc);
}

size_t Molecule::dimensionality(double eps) const {
    using chemfiles::Vector3D;
    auto& pos = frame_.positions();

    if (size() == 0 || size() == 1) {
        return 0;
    }

    // Any two points are on the same line
    if (size() == 2) {
        return 1;
    }

    Vector3D lin_vec = pos[1] - pos[0];
    Vector3D plane_vec = pos[2] - pos[0];

    if (size() == 3) {
        if (cos_sim(lin_vec, plane_vec) < eps) {
            return 1; // It's still linear!
        } else {
            return 2; // Any three points makes a plane
        }
    }

    Vector3D norm_vec = chemfiles::cross(lin_vec, plane_vec);
    auto d = chemfiles::dot(norm_vec, pos[0]);

    bool is_linear = true;
    for (size_t i = 3; i < size(); ++i) {
        Vector3D curr_vec = pos[i] - pos[0];
        if (is_linear && cos_sim(lin_vec, curr_vec) > eps) {
            is_linear = false;
        }

        if (!is_linear) {
            auto test_d = chemfiles::dot(norm_vec, pos[i]);
            if (d < eps && test_d < eps) {
                continue;
            }

            if (d < eps && test_d > eps) {
                return 3;
            }

            if (std::fabs(std::fmod(test_d, d)) > eps) {
                return 3;
            }
        }
    }

    if (is_linear) {
        return 1;
    } else {
        return 2;
    }
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
