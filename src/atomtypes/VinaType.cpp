#include "spear/atomtypes/VinaType.hpp"
#include "spear/Molecule.hpp"

#include <map>
#include <string>
#include <locale>

using namespace Spear;

enum xs_type {
    C_H = 0,
    C_P = 1,
    N_P = 2,
    N_D = 3,
    N_A = 4,
    N_DA = 5,
    O_P = 6, // Unused???
    O_D = 7,
    O_A = 8,
    O_DA = 9,
    S_P = 10,
    P_P = 11,
    F_H = 12,
    Cl_H = 13,
    Br_H = 14,
    I_H = 15,
    Metal_D = 16,
    SKIP = 17, // Skip for calculations
};

const char* const xs_name[]{
    "C_H",
    "C_P",
    "N_P",
    "N_D",
    "N_A",
    "N_DA",
    "O_P",
    "O_D",
    "O_A",
    "O_DA",
    "S_P",
    "P_P",
    "F_H",
    "Cl_H",
    "Br_H",
    "I_H",
    "Metal_D",
    "SKIP", // Skip for calculations
};

const double xs_vdw_radii[] = {
    1.9, // C_H
    1.9, // C_P
    1.8, // N_P
    1.8, // N_D
    1.8, // N_A
    1.8, // N_DA
    1.7, // O_P
    1.7, // O_D
    1.7, // O_A
    1.7, // O_DA
    2.0, // S_P
    2.1, // P_P
    1.5, // F_H
    1.8, // Cl_H
    2.0, // Br_H
    2.2, // I_H
    1.2, // Metal_D
    0.0, // Skip
};

#define AMINO_ACID_CORE {"OXT",O_A},{"CA",C_P},{"C",C_P},{"O",O_A},{"N",N_D}

const std::map<std::string, std::map<std::string,xs_type>> xs_standard_residues {
    {"ALA",{AMINO_ACID_CORE,{"CB",C_H}}},
    {"ARG",{AMINO_ACID_CORE,{"CB",C_H},{"CD",C_H},{"CG",C_H},{"CZ",C_P},{"NE",N_D},{"NH1",N_D},{"NH2",N_D}}},
    {"ASN",{AMINO_ACID_CORE,{"CB",C_H},{"CG",C_P},{"ND2",N_D},{"OD1",O_A}}},
    {"ASP",{AMINO_ACID_CORE,{"CB",C_H},{"CG",C_P},{"OD1",O_A},{"OD2",O_A}}},
    {"CYS",{AMINO_ACID_CORE,{"CB",C_P},{"SG",S_P}}},
    {"GLN",{AMINO_ACID_CORE,{"CB",C_H},{"CD",C_P},{"CG",C_H},{"NE2",N_D},{"OE1",O_A}}},
    {"GLU",{AMINO_ACID_CORE,{"CB",C_H},{"CD",C_P},{"CG",C_H},{"OE1",O_A},{"OE2",O_A}}},
    {"GLY",{AMINO_ACID_CORE}},
    {"HIS",{AMINO_ACID_CORE,{"CB",C_H},{"CD2",C_P},{"CE1",C_P},{"CG",C_P},{"ND1",N_A},{"NE2",N_A}}},
    {"HIE",{AMINO_ACID_CORE,{"CB",C_H},{"CD2",C_P},{"CE1",C_P},{"CG",C_P},{"ND1",N_A},{"NE2",N_D}}},
    {"HID",{AMINO_ACID_CORE,{"CB",C_H},{"CD2",C_P},{"CE1",C_P},{"CG",C_P},{"ND1",N_D},{"NE2",N_A}}},
    {"HIP",{AMINO_ACID_CORE,{"CB",C_H},{"CD2",C_P},{"CE1",C_P},{"CG",C_P},{"ND1",N_D},{"NE2",N_D}}},
    {"ILE",{AMINO_ACID_CORE,{"CB",C_H},{"CD1",C_H},{"CG1",C_H},{"CG2",C_H}}},
    {"LEU",{AMINO_ACID_CORE,{"CB",C_H},{"CD1",C_H},{"CD2",C_H},{"CG",C_H}}},
    {"LYS",{AMINO_ACID_CORE,{"CB",C_H},{"CD",C_H},{"CE",C_P},{"CG",C_H},{"NZ",N_D}}},
    {"MET",{AMINO_ACID_CORE,{"CB",C_H},{"CE",C_H},{"CG",C_P},{"SD",S_P}}},
    {"PHE",{AMINO_ACID_CORE,{"CB",C_H},{"CD1",C_H},{"CD2",C_H},{"CE1",C_H},{"CE2",C_H},{"CG",C_H},{"CZ",C_H}}},
    {"PRO",{AMINO_ACID_CORE,{"CB",C_H},{"CD",C_P},{"CG",C_H}}},
    {"SER",{AMINO_ACID_CORE,{"CB",C_P},{"OG",O_DA}}},
    {"THR",{AMINO_ACID_CORE,{"CB",C_P},{"CG2",C_H},{"OG1",O_DA}}},
    {"TRP",{AMINO_ACID_CORE,{"CB",C_H},{"CD1",C_P},{"CD2",C_H},{"CE2",C_P},{"CE3",C_H},{"CG",C_H},{"CH2",C_H},{"CZ2",C_H},{"CZ3",C_H},{"NE1",N_D}}},
    {"TYR",{AMINO_ACID_CORE,{"CB",C_H},{"CD1",C_H},{"CD2",C_H},{"CE1",C_H},{"CE2",C_H},{"CG",C_H},{"CZ",C_P},{"OH",O_DA}}},
    {"VAL",{AMINO_ACID_CORE,{"CB",C_H},{"CG1",C_H},{"CG2",C_H}}},
};

static bool is_in_strong_resonance(const AtomVertex& av) {
    for (auto bond : av.bonds()) {
        auto neighbor = (av == bond.source()) ? bond.target() : bond.source();
        if ((neighbor.atomic_number() == Element::O || neighbor.atomic_number() == Element::N) &&
            bond.order() == Bond::DOUBLE) {
            return true;
        }
    }
    return false;
}

// Easy case, just check if carbon is bound to a heteroatom
static xs_type get_c_xs_type(const AtomVertex& av) {
    for (auto neighbor : av.neighbors()) {
        if (neighbor.atomic_number() != Element::C && 
            neighbor.atomic_number() != Element::H) {
            return xs_type::C_P;
        }
    }

    return xs_type::C_H;
}

// Hard case, there's three bond counts to get, which in turn involve checks
static xs_type get_n_xs_type(const AtomVertex& av) {
    if (av.degree() == 0) { // Typically Ammonia, assume protonated
        return xs_type::N_DA;
    } else if (av.degree() == 1) { // Nitrile, monoamine, etc
        auto neighbor = av[0];
        auto bond = av.br().bond(av, neighbor);

        if (bond.order() == Bond::TRIPLE) {
            return xs_type::N_A; // It can only accept - Nitrile
        } else if (bond.order() == Bond::DOUBLE) {
            return xs_type::N_D; // Imine, or guanidine
        } else {
            return xs_type::N_D; // AutoDOCK labels amines as Donor only
        }
    } else { // A linker, amide, guanidine, amine of some sort, etc,
        size_t sum_of_orders = 0;
        size_t sum_of_hydrogens = av.explicit_hydrogens();
        bool is_withdraw = false;
        for (auto neighbor : av.neighbors()) {
            auto bond = av.br().bond(av, neighbor);

            switch (bond.order()) {
            case Bond::SINGLE:
                sum_of_orders += 1;
                break;
            case Bond::DOUBLE:
                sum_of_orders += 2;
                break;
            case Bond::TRIPLE:
                sum_of_orders += 3;
                break;
            // Withdrawn cases.
            case Bond::AMIDE:
            case Bond::AROMATIC:
                sum_of_orders += 1;
                is_withdraw = true;
                break;
            // Something bizzare - give up
            default:
                return xs_type::N_P;
                break;
            }

            is_withdraw = is_withdraw || is_in_strong_resonance(neighbor);
        }

        if (sum_of_hydrogens >= 2) {
            return xs_type::N_D;
        }

        // Check if the nitrogen is positively charged and contains H
        // or X#N-X or X=N=X
        if (av.degree() > 3 || (av.degree() == 2 && sum_of_orders > 3)) {
            return sum_of_hydrogens == 0 ? xs_type::N_P : xs_type::N_D;
        }

        // From here on down, we know the degree must be 2 or 3
        if (sum_of_hydrogens >= 1) {
            // Must be H donating, but is accepting? check withdraw
            return is_withdraw ? xs_type::N_D : xs_type::N_DA;
        }

        // Must be no explicit hydrogens. We now check for implicit Hs

        // If only single bonds and degree 2, there's an implicit H
        if (sum_of_orders == 2 && av.degree() == 2) {
            // Linking Amine(not withdraw - Accept) / Amide(withdraw - just donate)
            // or impure heterocycle (aromatic and donate only)
            return (is_withdraw || av.is_aromatic()) ? xs_type::N_D
                                                     : xs_type::N_DA;
        }

        // no implicit hydrogen, we must only accept
        // regardless of withdrawn because we are an SP2 nitro of some form
        if (sum_of_orders == 3 && av.degree() == 2) {
            return xs_type::N_A;
        }

        // Degree must be 3!
        if (sum_of_orders == 3 && !is_withdraw) {
            return xs_type::N_A; // Linking SP2 nitrogen, aromatic, etc
        }
        // Fall through to Polar nitrogen
    }
    return xs_type::N_P;
}

// Medium case, two bond cases - but no additional checks.
static xs_type get_o_xs_type(const AtomVertex& av, uint64_t mode) {

    if (av.degree() == 0) { // Typically water
        if ((mode & VinaType::DONAR_WATER) && (mode & VinaType::ACCEPTOR_WATER)) {
            return xs_type::O_DA;
        }
        
        if (mode & VinaType::DONAR_WATER) {
            return xs_type::O_D;
        }

        if (mode & VinaType::ACCEPTOR_WATER) {
            return xs_type::O_A;
        }

        return xs_type::SKIP;
    } else if (av.degree() == 1) { // Carbonyl or alcohol
        auto neighbor = av[0];
        auto bond = av.br().bond(av, neighbor);

        if (bond.order() == Bond::DOUBLE) {
            return xs_type::O_A; // Carbonyl
        } else if (is_in_strong_resonance(neighbor)) {
            return xs_type::O_A; // Carboxylic acid, phosphate, sulfate
        } else {
            return xs_type::O_DA; // Alcohol
        }
    } else if (av.degree() == 2) { // alcohol or ether
        for (auto neighbor : av.neighbors()) {
            if (neighbor.atomic_number() == Element::H) {
                return xs_type::O_DA;
            }
        }
        return xs_type::O_A; // ether, phosphodiester, etc
    }
    // Oddity, not possible in RCSB
    return xs_type::O_P;
}

static xs_type get_xs_type(const AtomVertex& av, uint64_t mode) {

    auto residue = av.br().topology().residue_for_atom(av);
    if (residue) {
        // Allow explicit protonation to override defaults
        if (residue->name() == "HIS" && (av.name() == "ND1" || av.name() == "NE2")) {
            if (av.explicit_hydrogens() != 0) {
                return xs_type::N_D;
            }

            if (av.name() == "ND1" && (mode & VinaType::HIS_ND1_PROTONATED)) {
                return xs_type::N_D;
            }

            if (av.name() == "NE2" && (mode & VinaType::HIS_NE2_PROTONATED)) {
                return xs_type::N_D;
            }
        }

        auto sr = xs_standard_residues.find(residue->name());
        if (sr != xs_standard_residues.end()) {
            auto at = sr->second.find(av.name());
            if (at != sr->second.end()) {
                return at->second;
            }
        }
    }

    switch (av.atomic_number()) {
    // Handle 'easy' cases
    case Element::F:
        return xs_type::F_H;
    case Element::P:
        return xs_type::P_P;
    case Element::S:
        return xs_type::S_P;
    case Element::Cl:
        return xs_type::Cl_H;
    case Element::Br:
        return xs_type::Br_H;
    case Element::I:
        return xs_type::I_H;
    case Element::Mg:
    case Element::Ca:
    case Element::Mn:
    case Element::Fe:
    case Element::Zn:
        return xs_type::Metal_D;
    case Element::C:
        return get_c_xs_type(av);
    case Element::N:
        return get_n_xs_type(av);
    case Element::O:
        return get_o_xs_type(av, mode);
    default:
        return xs_type::SKIP;
    }
}

VinaType::VinaType(const Molecule& mol, uint64_t mode) : mol_(mol), mode_(mode) {
    reserve(mol.size());

    for (auto av : mol) {
        push_back(get_xs_type(av, mode_));
    }
}

size_t VinaType::add_atom(size_t new_idx) {
    push_back(get_xs_type(mol_[new_idx], mode_));
}

void VinaType::add_bond(size_t idx1, size_t idx2, Bond::Order /*unused*/) {
    (*this)[idx1] = get_xs_type(mol_[idx1], mode_);
    (*this)[idx2] = get_xs_type(mol_[idx2], mode_);
}

bool VinaType::is_hydrophobic(size_t xs) {
    return xs == xs_type::C_H  || xs == xs_type::F_H || xs == xs_type::Cl_H ||
           xs == xs_type::Br_H || xs == xs_type::I_H;
}

bool VinaType::is_acceptor(size_t xs) {
    return xs == xs_type::N_A || xs == xs_type::N_DA || xs == xs_type::O_A ||
           xs == xs_type::O_DA;
}

bool VinaType::is_donor(size_t xs) {
    return xs == xs_type::N_D  || xs == xs_type::N_DA || xs == xs_type::O_D ||
           xs == xs_type::O_DA || xs == xs_type::Metal_D;
}

bool VinaType::is_skip(size_t xs) {
    return xs == xs_type::SKIP;
}

template<> std::string Spear::atomtype_name_for_id<VinaType>(size_t id) {
    return xs_name[id];
}

template<> size_t Spear::atomtype_id_for_name<VinaType>(const std::string& name) {
    return Element::SymbolForName.at(name);
}

template<> size_t Spear::atomtype_id_count<VinaType>() {
    return xs_type::SKIP + 1;
}

template<> double Spear::van_der_waals<VinaType>(size_t id) {
    return xs_vdw_radii[id];
}
