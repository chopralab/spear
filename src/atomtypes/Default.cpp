#include "spear/atomtypes/Default.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

#include <map>
#include <string>
#include <locale>

using namespace Spear;

static void dec_hybridization(Hybridization& h) {
    switch(h) {
        case Hybridization::SP2: h = Hybridization::SP; break;
        case Hybridization::SP3: h = Hybridization::SP2; break;
    }
}

Default::Default(const Molecule& mol) : mol_(mol) {
    atom_types_.reserve(mol.size());
    hybridizations_.reserve(mol.size());

    for (auto av : mol) {
        auto atomic_number = av.atomic_number();

        atom_types_.push_back(atomic_number);
        switch(atomic_number) {
            case 6:
            case 32:
            case 7:
            case 8:
            case 9:
            case 15:
            case 16:
            case 17:
            case 35:
            case 53:
            break; //Only break, rest are continues
            // Hydrogens
            case 1:
            hybridizations_.push_back(Hybridization::FORCED);
            continue;
            // Boron and Aluminium
            case 5:
            case 13:
            hybridizations_.push_back(Hybridization::SP2);
            continue;
            // The rest
            default:
            hybridizations_.push_back(Hybridization::UNKNOWN);
            continue;
        }

        if (atomic_number == 15 || atomic_number == 16) {
            hybridizations_.push_back(static_cast<Hybridization>(av.neighbor_count()));
            continue;
        }

        hybridizations_.push_back(Hybridization::SP3);

        for (auto bond : av.bonds()) {

            if (bond.order() == chemfiles::Bond::DOUBLE) {
                dec_hybridization(hybridizations_.back());
            }

            // Just do it twice
            if (bond.order() == chemfiles::Bond::TRIPLE) {
                dec_hybridization(hybridizations_.back());
                dec_hybridization(hybridizations_.back());
            }
        }
    }
}

bool Default::is_aromatic(size_t atom_id) const {
    return false;
}

template<> std::string atomtype_name_for_id<Default>(size_t id) {
    return "NA";
}

template<> size_t atomtype_id_for_name<Default>(std::string name) {
    chemfiles::Atom a(name);
    return *a.atomic_number();
}

template<> size_t atomtype_id_count<Default>() {
    return 118;
}

template<> double van_der_waals<Default>(size_t id) {
    return 0.0;
}
