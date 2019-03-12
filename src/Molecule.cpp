#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/hawick_circuits.hpp>
#include <boost/graph/undirected_graph.hpp>

#include "spear/AtomType.hpp"
#include "spear/atomtypes/Default.hpp"

using Eigen::Vector3d;

using namespace Spear;

size_t AtomVertex::expected_bonds() const {
    auto atomic_num = atomic_number();

    switch(atomic_num) {
        case Element::H:  // Hydrogen
        case Element::F:  // Fluorine
        case Element::Cl: // Chlorine
        case Element::Br: // Bromine
        case Element::I: // Iodine
        return 1;
        break;

        case Element::B:  // Boron
        case Element::Al: // Aluminium
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
struct cycle_saver {

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

    for (auto atom : frame_) {
        auto atomic_number = atom.atomic_number();
        if (atomic_number) {
            auto an = static_cast<Element::Symbol>(*atomic_number);
            boost::add_vertex(an, graph_);
        } else {
            boost::add_vertex(Element::Symbol(0), graph_);
        }
    }

    for (size_t i = 0; i < frame_.size(); ++i) {
        auto& pos = frame_.positions()[i];
        positions_.push_back({pos[0], pos[1], pos[2]});
    }

    // Create the graph representation
    auto bonds = topo.bonds();
    for ( size_t i = 0; i < bonds.size(); ++i) {
        boost::add_edge(bonds[i][0], bonds[i][1],
            static_cast<Bond::Order>(topo.bond_orders()[i]), graph_);
    }

    add_atomtype<Default>();
}

const std::set<std::set<size_t>> Molecule::rings() const {
    std::set<std::set<size_t>> ret_rings;

    cycle_saver vis(ret_rings);
    boost::hawick_circuits(graph_, vis);

    return ret_rings;
}

static double cos_sim(const Eigen::Vector3d& u, const Eigen::Vector3d& v) {
    auto arc = u.dot(v) / (u.norm() * v.norm());
    return arc >= 1? 0 : std::acos(arc);
}

size_t Molecule::dimensionality(double eps) const {
    using Eigen::Vector3d;

    if (size() == 0 || size() == 1) {
        return 0;
    }

    // Any two points are on the same line
    if (size() == 2) {
        return 1;
    }

    Vector3d lin_vec = positions_[1] - positions_[0];
    Vector3d plane_vec = positions_[2] - positions_[0];

    if (size() == 3) {
        if (cos_sim(lin_vec, plane_vec) < eps) {
            return 1; // It's still linear!
        } else {
            return 2; // Any three points makes a plane
        }
    }

    Vector3d norm_vec =lin_vec.cross(plane_vec);
    auto d = norm_vec.dot(positions_[0]);

    bool is_linear = true;
    for (size_t i = 3; i < size(); ++i) {
        Vector3d curr_vec = positions_[i] - positions_[0];
        if (is_linear && cos_sim(lin_vec, curr_vec) > eps) {
            is_linear = false;
        }

        if (!is_linear) {
            auto test_d = norm_vec.dot(positions_[i]);
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

std::vector<Spear::BondEdge> Molecule::get_bonds_in(const std::set<size_t>& atoms) const {
    auto all_bonds = bonds();

    std::vector<BondEdge> ret;

    std::copy_if(all_bonds.begin(), all_bonds.end(),
                 std::back_inserter(ret),
        [&atoms, this](BondEdge e){
            auto target = e.target();
            auto source = e.source();
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

void Molecule::remove_hydrogens() {
    auto& mol = *this;

    // Removal of an atom invalidates the vertex descriptors as we are using a
    // vecS for vertex storage (required to link index to frame).
    // Therefore we need to restart the search when we remove something.
    bool has_hydrogens = false;
    do {
        auto verticies_iter = boost::vertices(graph_);    
        has_hydrogens = false;
        for (auto v = verticies_iter.first; v != verticies_iter.second; ++v) {
            auto index = boost::get(boost::vertex_index_t(), graph_, *v);
            if (mol[index].atomic_number() == 1) {
                boost::clear_vertex(*v, graph_);
                boost::remove_vertex(*v, graph_);
                frame_.remove(index);
                has_hydrogens = true;
                break;
            }
        }
    } while (has_hydrogens);
}

AtomVertex Molecule::add_atom(Element::Symbol n_atom, const Eigen::Vector3d& pos) {
    auto desc = boost::add_vertex(n_atom, graph_);
    positions_.push_back(pos);
    frame_.add_atom(chemfiles::Atom(Element::Name[n_atom]), {pos[0], pos[1], pos[2]});
    for (const auto& at : atom_types_) {
        at.second->add_atom(desc);
    }
    return {this, desc};
}

BondEdge Molecule::add_bond(size_t idx1, size_t idx2, Bond::Order order) {
    auto new_edge = boost::add_edge(idx1, idx2, order, graph_);
    frame_.add_bond(idx1, idx2, static_cast<chemfiles::Bond::BondOrder>(order));
    return {this, new_edge.first};
}

static Eigen::Vector3d perpendicular(const Eigen::Vector3d& other) {
    Vector3d res(0.0, 0.0, 0.0);
    if (other[0]) {
      if (other[1]) {
        res[1] = -1 * other[0];
        res[0] = other[1];
      } else if (other[2]) {
        res[2] = -1 * other[0];
        res[0] = other[2];
      } else {
        res[1] = 1;
      }
    } else if (other[1]) {
      if (other[2]) {
        res[2] = -1 * other[1];
        res[1] = other[2];
      } else {
        res[0] = 1;
      }
    } else if (other[2]) {
      res[0] = 1;
    }
    res.normalize();
    return res;
}

// This code is largely based off of RDkit's code used to generate hydrogen
// locations. The original version is availible here:
// https://github.com/rdkit/rdkit/blob/95eeec7b3424f482416aa0c85b77fc73ab3008c1/Code/GraphMol/AddHs.cpp#L45

AtomVertex Molecule::add_atom_to(Element::Symbol n_atom, size_t index) {

    auto av = (*this)[index];

    Vector3d dirVect(0, 0, 0);

    Vector3d heavyPos = av.position();

    Vector3d perpVect, rotnAxis, nbrPerp;
    Vector3d nbr1Vect, nbr2Vect, nbr3Vect;

    auto dimensional = dimensionality();

    // All the coordinates are the same, let's ignore this case
    if (dimensional == 0 && av.neighbor_count() >= 1) {
        auto h_atom = add_atom(n_atom, heavyPos);
        add_bond(index, h_atom);
        return h_atom;
    }

    auto hybrid = get_default_atomtype()->hybridization(index);

    if (hybrid == Hybridization::UNKNOWN || hybrid == Hybridization::FORCED) {
        throw std::runtime_error("Bad hybridization state of atom!");
    }

    bool is3D = dimensional == 3;

    double bondLength =
        is3D ?
            Element::CovalentRadius[n_atom] +
            Element::CovalentRadius[av.atomic_number()] :
            1.0;

    switch (av.neighbor_count()) {
    case 0: // No other atom neighbors, add position in the z-direction
        dirVect[2] = 1;
        break;

    case 1: // One Neighbor present
        heavyPos = av.position();

        // get a normalized vector pointing away from the neighbor:
        nbr1Vect = av[0].position() - heavyPos;

        nbr1Vect.normalize();
        nbr1Vect *= -1;

        // nbr1Vect points away from the other atom, add new H accordingly:
        switch (hybrid) {
        case Hybridization::SP3:
            // get a perpendicular to nbr1Vect:
            if (is3D) {
                perpVect = perpendicular(nbr1Vect);
            } else {
                perpVect[2] = 1.0;
            }
            // and move off it:
            dirVect = Eigen::AngleAxisd((180 - 109.471) * M_PI / 180.0, perpVect) * nbr1Vect;
            break;
        case Hybridization::SP2:
            // default position is to just take an arbitrary perpendicular:
            perpVect = perpendicular(nbr1Vect);
            if (av[0].neighbor_count() > 1) {
                // can we use the neighboring atom to establish a perpendicular?
                auto nbrBond = *(av.bonds().begin()); // Only one bond, so this has to be the one!
                if (nbrBond.order() == Bond::AROMATIC ||
                    nbrBond.order() == Bond::DOUBLE) {
                    for (auto nbr2 : av[0].neighbors()) {
                        if (nbr2 != av) {
                            nbr2Vect = nbr2.position() - av[0].position();
                            nbr2Vect.normalize();
                            break;
                        }
                    }
                    perpVect = nbr2Vect.cross(nbr1Vect);
                }
            }
            perpVect.normalize();
            // rotate the nbr1Vect 60 degrees about perpVect and we're done:
            dirVect = Eigen::AngleAxisd(60.0 * M_PI / 180.0, perpVect) * nbr1Vect;
            break;
        case Hybridization::SP:
            // just lay the H along the vector:
            dirVect = nbr1Vect;
            break;
        default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the vector:
            dirVect = nbr1Vect;
        }
        break;
    case 2: // Two other neighbors:
        // start along the average of the two vectors:
        heavyPos = av.position();
        nbr1Vect = heavyPos - av[0].position();
        nbr2Vect = heavyPos - av[1].position();

        nbr1Vect.normalize();
        nbr2Vect.normalize();
        dirVect = nbr1Vect + nbr2Vect;

        dirVect.normalize();
        if (is3D) {
            switch (hybrid) {
            case Hybridization::SP3:
                // get the perpendicular to the neighbors:
                nbrPerp = nbr1Vect.cross(nbr2Vect);
                // and the perpendicular to that:
                rotnAxis = nbrPerp.cross(dirVect);
                // and then rotate about that:
                rotnAxis.normalize();
                dirVect = Eigen::AngleAxisd((109.471 / 2) * M_PI / 180., rotnAxis) * dirVect;
                break;
            case Hybridization::SP2:
                // don't need to do anything here, the H atom goes right on the direction vector
                break;
            case Hybridization::SP:
                // Can't add an atom here, the bond is saturated!
                throw std::runtime_error("Cannot add a third atom to SP atoms!");
                break;
            default:
                // FIX: handle other hybridizations
                // for now, just lay the H along the neighbor vector;
                break;
            }
        }
        // don't need to do anything here, the H atom goes right on the direction vector
        break;
    case 3: // Three other neighbors:

        if (hybrid == Hybridization::SP2 || hybrid == Hybridization::SP) {
            throw std::runtime_error("Cannot add a fourth atom to SP or SP2 atoms!");
        }

        // use the average of the three vectors:
        heavyPos = av.position();
        nbr1Vect = heavyPos - av[0].position();
        nbr2Vect = heavyPos - av[1].position();
        nbr3Vect = heavyPos - av[2].position();
        nbr1Vect.normalize();
        nbr2Vect.normalize();
        nbr3Vect.normalize();

        // if three neighboring atoms are more or less planar, this
        // is going to be in a quasi-random (but almost definitely bad)
        // direction...
        if (is3D) {
            dirVect = nbr1Vect + nbr2Vect + nbr3Vect;
        } else {
            // we're in flatland
            // We're in a 2D conformation, put the H between the two neighbors
            // that have the widest angle between them:
            double minDot = nbr1Vect.dot(nbr2Vect);
            dirVect = nbr1Vect + nbr2Vect;
            if (nbr2Vect.dot(nbr3Vect) < minDot) {
                minDot = nbr2Vect.dot(nbr3Vect);
                dirVect = nbr2Vect + nbr3Vect;
            }
            if (nbr1Vect.dot(nbr3Vect) < minDot) {
                minDot = nbr1Vect.dot(nbr3Vect);
                dirVect = nbr1Vect + nbr3Vect;
            }
            dirVect *= -1;
        }
        break;
    default: // Give up!
        break;
    }

    auto h_atom = add_atom(n_atom, heavyPos + dirVect * bondLength);
    add_bond(index, h_atom);

    return h_atom;
}
