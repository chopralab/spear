#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/hawick_circuits.hpp>
#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/connected_components.hpp>

#include "spear/AtomType.hpp"
#include "spear/atomtypes/Default.hpp"

#include "spear/Geometry.hpp"

#ifndef M_PI
static const auto M_PI = std::acos(0.0) * 2;
#endif

using Eigen::Vector3d;

using namespace Spear;

size_t AtomVertex::expected_bonds() const {
    auto atomic_num = atomic_number();

    switch(atomic_num) {
        case Element::H:  // Hydrogen
        case Element::F:  // Fluorine
        case Element::Cl: // Chlorine
        case Element::Br: // Bromine
        case Element::I:  // Iodine
        return 1;
        break;

        case Element::B:  // Boron
        case Element::Al: // Aluminium
        return 3;
        break;

        default: break;
    }

    auto hybridization_state = br_->get_default_atomtype()->hybridization(index_);

    if(atomic_num == Element::C) {
        switch (hybridization_state) {
            case Hybridization::SP: return 2; break;
            case Hybridization::SP2: return 3; break;
            case Hybridization::SP3: return 4; break;
            default: return 0; break;
        }
    }

    if(atomic_num == Element::N) {
        switch (hybridization_state) {
            case Hybridization::SP: return 1; break;
            case Hybridization::SP2: return 2; break;
            case Hybridization::SP3: return 3; break;
            default: return 0; break;
        }
    }

    if(atomic_num == Element::O) {
        switch (hybridization_state) {
            case Hybridization::SP: return 1; break;
            case Hybridization::SP2: return 1; break;
            case Hybridization::SP3: return 2; break;
            default: return 0; break;
        }
    }

    if(atomic_num == 15 || atomic_num == 16) { // Phosphorus and Sulfur
        if (is_aromatic()) return 2;
        return degree(); // Hard to tell, but things shouldn't be protonated
    }

    // Give up....
    return static_cast<size_t>(hybridization_state);
}

// The cycle_printer is a visitor that will print the path that comprises
// the cycle. Note that the back() vertex of the path is not the same as
// the front(). It is implicit in the listing of vertices that the back()
// vertex is connected to the front().
struct cycle_saver {

    cycle_saver(RingSet& counter):
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

    RingSet& found_rings;
};

void Molecule::init_(const chemfiles::Frame& frame) {
    for (const auto& atom : frame) {
        auto atomic_number = atom.atomic_number();
        if (atomic_number) {
            auto an = static_cast<Element::Symbol>(*atomic_number);
            boost::add_vertex(an, graph_);
        } else {
            boost::add_vertex(Element::Symbol(0), graph_);
        }
    }

    for (size_t i = 0; i < frame.size(); ++i) {
        auto& pos = frame.positions()[i];
        positions_.emplace_back(pos[0], pos[1], pos[2]);
    }

    // Create the graph representation
    auto bonds = topology_.bonds();
    auto bos = topology_.bond_orders();
    for ( size_t i = 0; i < bonds.size(); ++i) {
        boost::add_edge(bonds[i][0], bonds[i][1],
            static_cast<Bond::Order>(bos[i]), graph_);
    }

    add_atomtype<Default>();
}

void Molecule::rings_() {

    if (!all_rings_.empty()) {
        return;
    }

    cycle_saver vis(all_rings_);
    boost::hawick_circuits(graph_, vis);

    for (auto& ring : all_rings_) {
        for (auto& atom : ring) {
            atom_to_ring_.insert({atom, ring});
        }
    }
}

void Molecule::smallest_set_of_smallest_rings_() {

    if (!sssr_.empty()) {
        return;
    }

    auto all_rings = rings();

    if (all_rings.size() == 0) {
        return;
    }

    std::vector<int> component(size());
    auto num = static_cast<size_t>(boost::connected_components(graph_, &component[0]));
    auto sssr_count = boost::num_edges(graph_) + 1 + (num - 1) - size();

    if (sssr_count >= all_rings.size()) {
        sssr_ = all_rings_;
        return;
    }

    size_t max_degree = 0;
    for (auto av : *this) {
        max_degree = std::max(max_degree, av.degree());
    }

    auto find_independant_rings_with_degree =
    [&all_rings, this] (size_t degree) {
        std::vector<bool> used_nodes(size(), false); // move out for speed?

        // The rings are already sorted by their size, therefor we will
        // traverse all of the rings from the smallest to largest here.
        for (auto& current_ring : all_rings) {
            bool is_independant = false;

            // Have we seen this node before?
            // ... but only check if the current degree is correct!
            for (auto atom_id_search : current_ring) {
                if ((*this)[atom_id_search].degree() != degree ||
                    used_nodes[atom_id_search]) {
                    continue;
                }

                is_independant = true;
                break;
            }

            // This ring is completly encircled for the current degree
            if (!is_independant) {
                continue;
            }

            // We found a smallest ring! Mark its atoms as used
            sssr_.insert(current_ring);

            for (auto atom_id : current_ring) {
                used_nodes[atom_id] = true;
            }
        }
    };

    size_t iterations = 0; // I'm really scared of non termination. Let's be safe
    do {
        for (size_t i = 2; i <= max_degree; ++i) {
            find_independant_rings_with_degree(i);

            if (sssr_.size() == sssr_count) { // we're done!
                for (auto& ring : sssr_) {
                    for (auto& atom : ring) {
                        atom_to_sssr_.insert({atom, ring});
                    }
                }

                return;
            }
        }

        if (sssr_.size() > sssr_count) {
            break;
        }

        // things just got really weird - suppress the current SSSRs and try again
        for (const auto& added : sssr_) {
            all_rings.erase(added);
        }

        ++iterations;
    } while(all_rings.size() != 0 && iterations < 100);

    throw std::runtime_error(std::string("Unable to find SSSR: found ") +
                             std::to_string(sssr_.size()) + " rings, but expected " +
                             std::to_string(sssr_count));
}

static double cos_diff(const Eigen::Vector3d& u, const Eigen::Vector3d& v) {
    return 1.0 - std::abs(u.dot(v) / (u.norm() * v.norm()));
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

    if (size() == 3) {
        Vector3d lin_vec = positions_[1] - positions_[0];
        Vector3d plane_vec = positions_[2] - positions_[0];
        if (cos_diff(lin_vec, plane_vec) < eps) {
            return 1; // It's still linear!
        } else {
            return 2; // Any three points makes a plane
        }
    }

    // These two points, along with the first point, must be nonlinear
    // this allows one to calculate the 
    const Vector3d* nonlin1;
    const Vector3d* nonlin2;

    bool is_linear = true;
    for (size_t i = 3; i < size(); ++i) {
        Vector3d curr_vec = positions_[i] - positions_[0];
        if (is_linear) {
            // Check all previous points for non-linearity
            // If we find something non-linear, we store this point
            // and the other non-linear point to check for non-planarity
            for (size_t j = 0; j < i; ++j) {
                if (cos_diff(positions_[j] - positions_[0], curr_vec) > eps) {
                    nonlin1 = &positions_[j];
                    nonlin2 = &positions_[i];
                    is_linear = false;
                    break;
                }
            }

            continue;
        }

        auto nplane = nonplanar(positions_[0], positions_[i], *nonlin1, *nonlin2);
        if (std::abs(nplane) > eps) {
            return 3;
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

void Molecule::remove_bond(size_t idx1, size_t idx2) {
    boost::remove_edge(idx1, idx2, graph_);
    topology_.remove_bond(idx1, idx2);
}

void Molecule::remove_hydrogens() {
    auto& mol = *this;

    // Removal of an atom invalidates the vertex descriptors as we are using a
    // vecS for vertex storage (required to link index to frame).
    // Therefore we need to restart the search when we remove something.
    bool has_hydrogens = false;
    size_t skip_len = 0;
    do {
        auto verticies_iter = boost::vertices(graph_);    
        has_hydrogens = false;
        for (auto v = verticies_iter.first + static_cast<std::ptrdiff_t>(skip_len);
                  v != verticies_iter.second; ++v) {
            auto index = boost::get(boost::vertex_index, graph_, *v);
            if (mol[index].atomic_number() == Element::H) {
                skip_len = *v;
                boost::clear_vertex(*v, graph_);
                boost::remove_vertex(*v, graph_);
                topology_.remove(index);
                positions_.erase(positions_.begin() + index);
                for (const auto& at : atom_types_) {
                    at.second->remove_atom(index);
                }
                has_hydrogens = true;
                break;
            }
        }
    } while (has_hydrogens);
}

AtomVertex Molecule::add_atom(Element::Symbol n_atom, const Eigen::Vector3d& pos) {
    auto desc = boost::add_vertex(n_atom, graph_);
    positions_.push_back(pos);
    topology_.add_atom(chemfiles::Atom(Element::Name[n_atom]));
    for (const auto& at : atom_types_) {
        at.second->add_atom(desc);
    }
    return {this, desc};
}

BondEdge Molecule::add_bond(size_t idx1, size_t idx2, Bond::Order order) {
    auto new_edge = boost::add_edge(idx1, idx2, order, graph_);
    if (!new_edge.second) {
        throw std::runtime_error("Could not add bond between: " +
                                 std::to_string(idx1) + " and " +
                                 std::to_string(idx2));
    }
    topology_.add_bond(idx1, idx2, static_cast<chemfiles::Bond::BondOrder>(order));
    for (const auto& at : atom_types_) {
        at.second->add_bond(idx1, idx2, order);
    }
    return {this, new_edge.first};
}

static Eigen::Vector3d perpendicular(const Eigen::Vector3d& other) {
    Vector3d res(0.0, 0.0, 0.0);
    if (other[0] != 0.0) {
        if (other[1] != 0.0) {
            res[1] = -1 * other[0];
            res[0] = other[1];
        } else if (other[2] != 0.0) {
            res[2] = -1 * other[0];
            res[0] = other[2];
        } else {
            res[1] = 1;
        }
    } else if (other[1] != 0.0) {
        if (other[2] != 0.0) {
            res[2] = -1 * other[1];
            res[1] = other[2];
        } else {
            res[0] = 1;
        }
    } else if (other[2] != 0.0) {
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
    if (dimensional == 0 && av.degree() >= 1) {
        auto h_atom = add_atom(n_atom, heavyPos);
        add_bond(index, h_atom);
        return h_atom;
    }

    auto hybrid = get_default_atomtype()->hybridization(index);

    if (hybrid == Hybridization::UNKNOWN || hybrid == Hybridization::FORCED) {
        throw std::runtime_error("Cannot add atom to UNKNOWN or FORCED hybridizations");
    }

    // Spoof SP2 for planar nitrogens so that atoms are added in the same plane
    // and we can prevent the loss of conjugation
    if (av.atomic_number() == Element::N && hybrid == Hybridization::SP3 &&
        get_default_atomtype()->is_planar(index)) {
        hybrid = Hybridization::SP2;
    }

    bool is3D = dimensional == 3;

    double bondLength =
        is3D ?
            Element::CovalentRadius[n_atom] +
            Element::CovalentRadius[av.atomic_number()] :
            1.0;

    if (is3D && av.atomic_number() == Element::C && hybrid == Hybridization::SP2) {
        bondLength += 0.04;
    }

    if (is3D && av.atomic_number() == Element::C && hybrid == Hybridization::SP3) {
        bondLength += 0.07;
    }

    switch (av.degree()) {
    case 0: // No other atom neighbors, add position in the z-direction
        dirVect[2] = 1;
        break;

    case 1: // One Neighbor present
        // get a normalized vector pointing away from the neighbor:
        nbr1Vect = av[0].position() - heavyPos;

        nbr1Vect.normalize();
        nbr1Vect *= -1;

        // nbr1Vect points away from the other atom, add new H accordingly:
        switch (hybrid) {
        case Hybridization::SP3:
            // get a perpendicular to nbr1Vect:
            perpVect = perpendicular(nbr1Vect);
            if (!is3D) {
                perpVect[0] = 0.0;
                perpVect[1] = 0.0;
                perpVect[2] = 1.0;
            }
            // and move off it:
            dirVect = Eigen::AngleAxisd((180 - 109.471) * M_PI / 180.0, perpVect) * nbr1Vect;
            break;
        case Hybridization::SP2:
            // default position is to just take an arbitrary perpendicular:
            perpVect = perpendicular(nbr1Vect);
            if (av[0].degree() > 1) {
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
            throw std::runtime_error("Cannot add a fourth atom to planar nitrogens or SP/SP2 atoms!");
        }

        // use the average of the three vectors:
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
            dirVect.normalize();
        } else {
            // we're in flatland
            // We're in a 2D conformation, put the H between the two neighbors
            // that have the widest angle between them:
            double minDot = nbr1Vect.dot(nbr2Vect);
            dirVect = nbr1Vect + nbr2Vect;
            if (nbr2Vect.dot(nbr3Vect) < minDot) {
                nbr2Vect.dot(nbr3Vect);
                dirVect = nbr2Vect + nbr3Vect;
            }
            if (nbr1Vect.dot(nbr3Vect) < minDot) {
                nbr1Vect.dot(nbr3Vect);
                dirVect = nbr1Vect + nbr3Vect;
            }
            dirVect *= -1;
        }
        break;
    case 4:
        if (hybrid == Hybridization::SP3 ||
            hybrid == Hybridization::SP2 ||
            hybrid == Hybridization::SP) {
            throw std::runtime_error("Cannot add a fifth atom to SP, SP2, or SP3 atoms!");
        }

    default: // Give up!
        break;
    }

    dirVect.normalize();
    auto h_atom = add_atom(n_atom, heavyPos + dirVect * bondLength);
    add_bond(av, h_atom);

    const auto& res = topology_.residue_for_atom(av);
    if (res) {
        chemfiles::Residue res2(res->name());
        res2.add_atom(h_atom);
        auto neigh_name = av.name();
        neigh_name[0] = 'H';
        topology_.add_residue(res2);
        topology_[h_atom].set_name(neigh_name);
    }

    return h_atom;
}

size_t Molecule::add_hydrogens() {
    bool hydrogens_to_add = false;
    size_t hydrogens_added = 0;
    do {
        hydrogens_to_add = false;
        for (auto atom : *this) {
            if (atom.implicit_hydrogens() < 1 || atom.atomic_number() == Element::H) {
                continue;
            }
            hydrogens_to_add = true;
            ++hydrogens_added;
            add_atom_to(Element::H, atom);
            break;
        }
    } while (hydrogens_to_add);

    return hydrogens_added;
}
