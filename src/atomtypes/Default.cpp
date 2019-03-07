#include "spear/atomtypes/Default.hpp"
#include "spear/Molecule.hpp"
#include <map>
#include <string>
#include <locale>
#include <iostream>
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
            case 1:
            hybridizations_.push_back(Hybridization::FORCED);
            break;
            case 15:
            case 16:
            hybridizations_.push_back(static_cast<Hybridization>(av.neighbor_count()));
            break;
            case 5:
            case 13:
            hybridizations_.push_back(Hybridization::SP2);
            break;
            case 6:
            case 32:
            case 7:
            case 8:
            case 9:
            case 17:
            case 35:
            case 53:
            hybridizations_.push_back(Hybridization::SP3);
            break;
            default:
            hybridizations_.push_back(Hybridization::UNKNOWN);
            break;
        }
    }

    auto& bonds = mol_.frame().topology().bonds();
    auto& bond_orders = mol_.frame().topology().bond_orders();

    for (size_t i = 0; i < bonds.size(); ++i) {
        if (bond_orders[i] == chemfiles::Bond::DOUBLE) {
            if (mol_[bonds[i][0]].atomic_number() != 15 &&
                mol_[bonds[i][0]].atomic_number() != 16) {
                std::cout << bonds[i][0] << mol_[bonds[i][0]].atomic_number() << " " << hybridizations_[bonds[i][0]] << " ";
                dec_hybridization(hybridizations_[bonds[i][0]]);
                std::cout << hybridizations_[bonds[i][0]] << std::endl;
            }
            if (mol_[bonds[i][1]].atomic_number() != 15 &&
                mol_[bonds[i][1]].atomic_number() != 16) {
                dec_hybridization(hybridizations_[bonds[i][1]]);
            }
        }

        // Just do it twice
        if (bond_orders[i] == chemfiles::Bond::TRIPLE) {
            if (mol_[bonds[i][0]].atomic_number() != 15 &&
                mol_[bonds[i][0]].atomic_number() != 16) {
                dec_hybridization(hybridizations_[bonds[i][0]]);
                dec_hybridization(hybridizations_[bonds[i][0]]);
            }
            if (mol_[bonds[i][1]].atomic_number() != 15 &&
                mol_[bonds[i][1]].atomic_number() != 16) {
                dec_hybridization(hybridizations_[bonds[i][1]]);
                dec_hybridization(hybridizations_[bonds[i][1]]);
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
