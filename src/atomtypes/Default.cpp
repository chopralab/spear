#include "spear/atomtypes/Default.hpp"
#include "spear/Molecule.hpp"

#include <map>
#include <string>
#include <locale>

using namespace Spear;

static void dec_hybridization(Hybridization& h) {
    switch(h) {
        case Hybridization::SP2: h = Hybridization::SP; break;
        case Hybridization::SP3: h = Hybridization::SP2; break;
        default: break; // Do nothing
    }
}

Default::Default(const Molecule& mol) : mol_(mol) {
    reserve(mol.size());
    hybridizations_.reserve(mol.size());

    aromatics_.resize(mol.size(), false);

    for (auto av : mol) {
        auto atomic_number = av.atomic_number();

        push_back(atomic_number);
        switch(atomic_number) {
        case Element::C:
        case Element::Si:
        case Element::N:
        case Element::O:
        case Element::F:
        case Element::P:
        case Element::S:
        case Element::Cl:
        case Element::Br:
        case Element::I:
            break; //Only break, rest are continues

        case Element::H:
            hybridizations_.push_back(Hybridization::FORCED);
            continue;

        case Element::B:
        case Element::Al:
            hybridizations_.push_back(Hybridization::SP2);
            continue;

        // The rest
        default:
            hybridizations_.push_back(Hybridization::UNKNOWN);
            continue;
        }

        if (atomic_number == 15 || atomic_number == 16) {
            hybridizations_.push_back(static_cast<Hybridization>(av.degree()));
            continue;
        }

        hybridizations_.push_back(Hybridization::SP3);

        bool is_aromatic = false; // used to prevent dec_hybridization twice

        for (auto bond : av.bonds()) {

            if (bond.order() == Bond::DOUBLE) {
                dec_hybridization(hybridizations_.back());
            }

            if (bond.order() == Bond::AROMATIC && !is_aromatic) {
                is_aromatic = true;
                dec_hybridization(hybridizations_.back());
            }

            // Just do it twice
            if (bond.order() == Bond::TRIPLE) {
                dec_hybridization(hybridizations_.back());
                dec_hybridization(hybridizations_.back());
            }
        }
    }

    assign_aromaticity();
    assign_aromaticity();
}

bool Default::is_aromatic(size_t atom_id) const {
    return aromatics_[atom_id];
}

static bool is_delocalized(const AtomVertex& av) {
    if (av.atomic_number()  == Element::S && av.degree() >= 3) {
        return false;
    }
    for (auto bond : av.bonds()) {
        switch(bond.order()) {
        case Bond::DOUBLE:
        case Bond::TRIPLE:
        case Bond::AROMATIC:
            return true;
        default:
            break;
        }
    }

    return false;
}

bool Default::is_planar(size_t atom_id) const {
    auto av = mol_[atom_id];

    if ( av.degree() <= 3 &&
        (av.atomic_number() == Element::B ||
         av.atomic_number() == Element::Al)) {
        return true;
    }

    for (auto bond : av.bonds()) {
        switch (bond.order()) {
        case Bond::AROMATIC:
        case Bond::TRIPLE: // Technically, a line is in a plane
            return true;
        case Bond::DOUBLE: // Here's where it gets tricky
            if (av.atomic_number() == Element::S &&
                av.degree() == 1) { // not a sulfoxide, sulfone, etc
                return true;
            } else if (av.atomic_number() == Element::P &&
                av.degree() <= 2) { // another odd case for phosphorus (eg P-oxide)
                return true;
            } else {
                return true;
            }
        case Bond::SINGLE: // Even worse, we need to check for delocalization for single LP groups
            if (av.atomic_number() == Element::N ||
                av.atomic_number() == Element::P) {
                if (av == bond.target() && is_delocalized(bond.source())) {
                    return true;
                }
                if (av == bond.source() && is_delocalized(bond.target())) {
                    return true;
                }
            }
            break;
        default:
            break;
        }
    }

    return false;
}

void Default::assign_aromaticity() {
    auto set_aromatics = [&] (const std::set<size_t>& r) {
        for (auto atom : r) {
            aromatics_[atom] = true;
        }
    };

    for (auto& sssr : mol_.smallest_set_of_smallest_rings()) {
        auto bonds_in = mol_.get_bonds_in(sssr);
        size_t num_explicit_aromatic = 0;
        size_t num_electrons = 0;
        for (const auto& b : bonds_in) {
            if (b.order() == Bond::AROMATIC) {
                num_explicit_aromatic++;
                num_electrons += 2;
            }
            if (b.order() == Bond::DOUBLE) {
                num_electrons += 2;
            }
        }

        // explicitly aromatic
        if (num_explicit_aromatic == sssr.size()) {
            set_aromatics(sssr);
            continue;
        }

        // Follow RDKit and don't assign aromaticity unless its a ring of 5/6
        if (sssr.size() != 5 && sssr.size() != 6) {
            continue;
        }

        bool canbe_aromatic = true;
        for (auto a : sssr) {
            switch (mol_[a].atomic_number()) {
            case Element::O:
            case Element::S:
            case Element::Te:
            case Element::Se:
                if (hybridizations_[a] == Hybridization::SP3 &&
                    mol_[a].degree() == 2) {
                    num_electrons += 2;
                } else {
                    canbe_aromatic = false;
                }
                break;
            case Element::N:
            case Element::P:
                // We've already accounted for double bonds
                if (hybridizations_[a] == Hybridization::SP3 &&
                    is_planar(a)) {
                    num_electrons += 2;
                }
                break;
            case Element::C:

                if (hybridizations_[a] != Hybridization::SP2) {
                    canbe_aromatic = false;
                }
                break;
            default:
                canbe_aromatic = false;
            }

            if (!canbe_aromatic) {
                break;
            }
        }

        // Found an element or condition that is not possibly aromatic
        if (!canbe_aromatic) {
            continue;
        }

        set_aromatics(sssr);
    }
}

size_t Default::add_atom(size_t new_idx) {
    auto av = mol_[new_idx];
    auto atomic_number = av.atomic_number();
    push_back(atomic_number);
    if (atomic_number == 1 || atomic_number == 2) {
        hybridizations_.push_back(Hybridization::FORCED);
        return atomic_number;
    }

    if (av.is_non_metal()) {
        hybridizations_.push_back(Hybridization::SP3);
        return atomic_number;
    }

    hybridizations_.push_back(Hybridization::UNKNOWN);
    return atomic_number;
}

void Default::add_bond(size_t idx1, size_t idx2, Bond::Order bo) {
    switch (bo) {
    case Bond::SINGLE:
    case Bond::AMIDE:
        return; // Do nothing
    case Bond::TRIPLE:
        dec_hybridization(hybridizations_[idx1]);
        dec_hybridization(hybridizations_[idx2]);
        // fall through
    case Bond::DOUBLE:
        dec_hybridization(hybridizations_[idx1]);
        dec_hybridization(hybridizations_[idx2]);
        break;
    case Bond::AROMATIC: // Benzyne?
        if ((*this)[idx1] == Element::C)
            hybridizations_[idx1] = Hybridization::SP2;
        if ((*this)[idx2] == Element::C)
            hybridizations_[idx2] = Hybridization::SP2;
        break;
    default:
        break;
    }
}

template<> std::string Spear::atomtype_name_for_id<Default>(size_t id) {
    return Element::Name[id];
}

template<> size_t Spear::atomtype_id_for_name<Default>(const std::string& name) {
    return Element::SymbolForName.at(name);
}

template<> size_t Spear::atomtype_id_count<Default>() {
    return 119; // Includes lone-pair
}

template<> double Spear::van_der_waals<Default>(size_t /*unused*/) {
    return 0.0;
}
