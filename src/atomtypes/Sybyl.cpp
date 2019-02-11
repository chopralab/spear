#include "spear/atomtypes/Sybyl.hpp"

#include "spear/Molecule.hpp"

using namespace Spear;

enum sybyl {
    C_3,   C_2, C_1, C_ar, C_cat, N_3, N_2, N_1, N_ar, N_am,
    N_pl3, N_4, O_3, O_2,  O_co2, S_3, S_2, S_o, S_o2, P_3,
    Ti_th, Ti_oh, Cr_th, Cr_oh, Co_oh, Ru_oh, H,

    Ac, Ag, Al, Am, Ar, As, At, Au, B,  Ba, Be, Bh, Bi, Bk, Br,
    Ca, Cd, Ce, Cf, Cm, Cs, Cu, Db, Ds, Dy, Er, Es, Eu, F,  Fe,
    Fm, Fr, Ga, Gd, Ge, He, Hf, Hg, Ho, Hs, I,  In, Ir, K,  Kr,
    La, Li, Lr, Lu, Lw, Md, Mg, Mn, Mo, Mt, Na, Nb, Nd, Ne, Ni,
    No, Np, Os, Pa, Pb, Pd, Pm, Po, Pr, Pt, Pu, Ra, Rb, Re, Rf,
    Rh, Rn, Sb, Sc, Se, Sg, Si, Sm, Sn, Sr, Ta, Tb, Tc, Te, Th,
    Tl, Tm, U,V, W, Xe, Y, Yb, Zn, Zr, Du
};

const std::unordered_map<std::string, size_t> sybyl_mask{
    {"C.3", C_3}, {"C.2", C_2}, {"C.1", C_1}, {"C.ar", C_ar}, {"C.cat", C_cat},
    {"N.3", N_3}, {"N.2", N_2}, {"N.1", N_1}, {"N.ar", N_ar}, {"N.am", N_am},
    {"N.pl3",N_pl3}, {"N.4", N_4}, {"O.3", O_3}, {"O.2", O_2},{"O.co2", O_co2},
    {"S.3", S_3}, {"S.2", S_2}, {"S.o", S_o}, {"S.o2", S_o2}, {"P.3", P_3},
    {"Ti.th",Ti_th}, {"Ti.oh",Ti_oh}, {"Cr.th", Cr_th}, {"Cr.oh",Cr_oh}, 
    {"Co.oh",Co_oh}, {"Ru.oh",Ru_oh}, {"H",H},

    {"Ac", Ac},   {"Ag", Ag},   {"Al", Al},   {"Am", Am},  {"Ar", Ar},
    {"As", As},   {"At", At},   {"Au", Au},   {"B", B},    {"Ba", Ba},
    {"Be", Be},   {"Bh", Bh},   {"Bi", Bi},   {"Bk", Bk},  {"Br", Br},
    {"Ca", Ca},   {"Cd", Cd},   {"Ce", Ce},   {"Cf", Cf},  {"Cm", Cm},
    {"Cs", Cs},   {"Cu", Cu},   {"Db", Db},
    {"Ds", Ds},   {"Dy", Dy},   {"Er", Er},   {"Es", Es},  {"Eu", Eu},
    {"F", F},     {"Fe", Fe},   {"Fm", Fm},   {"Fr", Fr},  {"Ga", Ga},
    {"Gd", Gd},   {"Ge", Ge},   {"He", He},
    {"Hf", Hf},   {"Hg", Hg},   {"Ho", Ho},   {"Hs", Hs},  {"I", I},
    {"In", In},   {"Ir", Ir},   {"K", K},     {"Kr", Kr},  {"La", La},
    {"Li", Li},   {"Lr", Lr},   {"Lu", Lu},   {"Lw", Lw},  {"Md", Md},
    {"Mg", Mg},   {"Mn", Mn},   {"Mo", Mo},   {"Mt", Mt},  {"Na", Na},
    {"Nb", Nb},   {"Nd", Nd},   {"Ne", Ne},   {"Ni", Ni},  {"No", No},
    {"Np", Np},   {"Os", Os},   {"Pa", Pa},   {"Pb", Pb},  {"Pd", Pd},
    {"Pm", Pm},   {"Po", Po},   {"Pr", Pr},   {"Pt", Pt},  {"Pu", Pu},
    {"Ra", Ra},   {"Rb", Rb},   {"Re", Re},   {"Rf", Rf},  {"Rh", Rh},
    {"Rn", Rn},   {"Sb", Sb},   {"Sc", Sc},   {"Se", Se},  {"Sg", Sg},
    {"Si", Si},   {"Sm", Sm},   {"Sn", Sn},   {"Sr", Sr},  {"Ta", Ta},
    {"Tb", Tb},   {"Tc", Tc},   {"Te", Te},   {"Th", Th},  {"Tl", Tl},
    {"Tm", Tm},   {"U", U},     {"V", V},     {"W", W},    {"Xe", Xe},
    {"Y", Y},     {"Yb", Yb},   {"Zn", Zn},   {"Zr", Zr},  {"Du", Du},
};

const char* const sybyl_unmask[] {
    "C.3",   "C.2", "C.1", "C.ar", "C.cat", "N.3", "N.2", "N.1", "N.ar", "N.am",
    "N.pl3", "N.4", "O.3", "O.2",  "O.co2", "S.3", "S.2", "S.o", "S.o2", "P.3",
    "Ti.th", "Ti.oh", "Cr.th", "Cr.oh", "Co.oh", "Ru.oh", "H",

    "Ac",  "Ag",  "Al",  "Am",  "Ar",  "As",  "At",  "Au",  "B",   "Ba",  "Be",
    "Bh",  "Bi",  "Bk",  "Br",  "Ca",  "Cd",  "Ce",  "Cf",  "Cm",  "Cs",  "Cu",
    "Db",  "Ds",  "Dy",  "Er",  "Es",  "Eu",  "F",   "Fe",  "Fm",  "Fr",  "Ga",
    "Gd",  "Ge",  "He",  "Hf",  "Hg",  "Ho",  "Hs",  "I",   "In",  "Ir",  "K",
    "Kr",  "La",  "Li",  "Lr",  "Lu",  "Lw",  "Md",  "Mg",  "Mn",  "Mo",  "Mt",
    "Na",  "Nb",  "Nd",  "Ne",  "Ni",  "No",  "Np",  "Os",  "Pa",  "Pb",  "Pd",
    "Pm",  "Po",  "Pr",  "Pt",  "Pu",  "Ra",  "Rb",  "Re",  "Rf",  "Rh",  "Rn",
    "Sb",  "Sc",  "Se",  "Sg",  "Si",  "Sm",  "Sn",  "Sr",  "Ta",  "Tb",  "Tc",
    "Te",  "Th",  "Tl",  "Tm",  "U",   "V",   "W",   "Xe",  "Y",   "Yb",  "Zn",
    "Zr",  "Du",
};

static size_t num_nonmetal(const Molecule& mol, AtomVertex atom) {
    size_t num_nonmetal = 0;
    for (auto a : atom) {
        auto element = mol[a].atomic_number();
        if (element <= 2 ||
            (element >= 5  && element <= 10) ||
            (element >= 14 && element <= 18) ||
            (element >= 32 && element <= 36) ||
            (element >= 51 && element <= 54) ){
            ++num_nonmetal;
        }
    }

    return num_nonmetal;
}

static size_t freeOxygens(const Molecule& mol, const Spear::AtomVertex& atom) {
    size_t freeOxygens = 0;
    for (auto neighbor : atom) {
        auto bondee = mol[neighbor];
        if (bondee.atomic_number() == 8 && num_nonmetal(mol, bondee) == 1) {
            ++freeOxygens;
        }
    }

    return freeOxygens;
}

static bool is_decloc(const Spear::AtomVertex& atom,
                      const std::vector<chemfiles::Bond>& bonds,
                      const std::vector<chemfiles::Bond::BondOrder>& bos) {
    for (size_t i = 0; i < bonds.size(); ++i) {
        if (bonds[i][0] != atom && bonds[i][1] != atom) {
            continue;
        }

        switch(bos[i]) {
            case chemfiles::Bond::DOUBLE:
            case chemfiles::Bond::TRIPLE:
            case chemfiles::Bond::AROMATIC:
                return true;
            default:
                break;
        }
    }

    return false;
}

Sybyl::Sybyl(const Molecule& mol, TypingMode mode) :
    atom_types_(mol.size(), sybyl::Du), mol_(mol) {
    switch (mode) {
    case GEOMETRY:
        name_ = "Sybyl_geometry";
        type_atoms_3d_();
        break;
    case TOPOLOGY:
        name_ = "Sybyl_topology";
        type_atoms_topo_();
        break;
    default:
        throw std::invalid_argument("Unknown mode");
        break;
    }
}

void Sybyl::type_atoms_3d_() {
    type_atoms_topo_();
    for (auto atom : mol_) {
        switch (atom.atomic_number()) {
        case 6:
            atom_types_[atom] = assign_carbon_3d_(atom);
            break;
        case 7:
            atom_types_[atom] = assign_nitrogen_3d_(atom);
            break;
        default:
            break;
        }
    }
}

void Sybyl::type_atoms_topo_() {
    for (auto atom : mol_) {
        switch (atom.atomic_number()) {
        case 1:
            atom_types_[atom] = H;
            break;
        case 6:
            atom_types_[atom] = assign_carbon_topo_(atom);
            break;
        case 7:
            atom_types_[atom] = assign_nitrogen_topo_(atom);
            break;
        case 8:
            atom_types_[atom] = assign_oxygen_(atom);
            break;
        case 16:
            atom_types_[atom] = assign_sulfur_(atom);
            break;
        case 15:
            atom_types_[atom] = P_3;
            break;
        case 27:
            atom_types_[atom] = Co_oh;
            break;
        case 44:
            atom_types_[atom] = Ru_oh;
            break;
        case 22:
            atom_types_[atom] = atom.neighbor_count() <= 4 ? Ti_th
                                                           : Ti_oh;
            break;
        case 24:
            atom_types_[atom] = atom.neighbor_count() <= 4 ? Cr_th
                                                           : Cr_oh;
            break;
        default:
            atom_types_[atom] = sybyl_mask.at(atom.type());
            break;
        }
    }
}

size_t Sybyl::assign_carbon_topo_(AtomVertex& atom){
    size_t num_double = 0, num_triple = 0, num_aromatic = 0;
    size_t num_nitrogen = 0;

    const auto& bonds = mol_.frame().topology().bonds();
    const auto& bond_orders = mol_.frame().topology().bond_orders();

    for (size_t i = 0; i < bonds.size(); ++i) {
        if (bonds[i][0] != atom && bonds[i][1] != atom) {
            continue;
        }

        if (mol_[bonds[i][0]].atomic_number() == 7) {
            num_nitrogen += (freeOxygens(mol_, mol_[bonds[i][0]]) == 0);
        }

        if (mol_[bonds[i][1]].atomic_number() == 7) {
            num_nitrogen += (freeOxygens(mol_, mol_[bonds[i][1]]) == 0);
        }

        if (bond_orders[i] == chemfiles::Bond::DOUBLE) {
            ++num_double;
        } else if (bond_orders[i] == chemfiles::Bond::TRIPLE) {
            ++num_triple;
        } else if (bond_orders[i] == chemfiles::Bond::AROMATIC) {
            ++num_aromatic;
        }
    }

    auto neighbor_count = atom.neighbor_count();

    if (neighbor_count >= 4 && num_double == 0 && num_triple == 0) {
        return C_3;
    }

    if (num_aromatic >= 2) {
        return C_ar;
    }

    if (num_nitrogen == 3 && neighbor_count == 3) {
        // Need to check for acyclic...
        return C_cat;
    }

    if (num_triple == 1 || num_double == 2) {
        return C_1;
    }

    return C_2;
}

size_t Sybyl::assign_nitrogen_topo_(AtomVertex& atom) {
    size_t num_double = 0, num_triple = 0, num_aromatic = 0;
    size_t num_amide = 0, num_deloc = 0;

    const auto& bonds = mol_.frame().topology().bonds();
    const auto& bond_orders = mol_.frame().topology().bond_orders();

    for (size_t i = 0; i < bonds.size(); ++i) {
        if (bonds[i][0] != atom && bonds[i][1] != atom) {
            continue;
        }

        if (mol_[bonds[i][0]].atomic_number() == 6) {
            num_amide += (freeOxygens(mol_, mol_[bonds[i][0]]) != 0);
            num_deloc += is_decloc(mol_[bonds[i][0]], bonds, bond_orders);
            size_t num_nitrogens = 0;
            for (auto bondee_bondee : mol_[bonds[i][0]]) {
                num_nitrogens += mol_[bondee_bondee].atomic_number() == 7;
            }

            if (num_nitrogens >= 3 && bond_orders[i] == 2) {
                return N_pl3;
            }
        }

        if (mol_[bonds[i][1]].atomic_number() == 6) {
            num_amide += (freeOxygens(mol_, mol_[bonds[i][1]]) != 0);
            num_deloc += is_decloc(mol_[bonds[i][1]], bonds, bond_orders);
            size_t num_nitrogens = 0;
            for (auto bondee_bondee : mol_[bonds[i][1]]) {
                num_nitrogens += mol_[bondee_bondee].atomic_number() == 7;
            }

            if (num_nitrogens >= 3 && bond_orders[i] == 2) {
                return N_pl3;
            }
        }

        if (bond_orders[i] == chemfiles::Bond::DOUBLE) {
            ++num_double;
        } else if (bond_orders[i] == chemfiles::Bond::TRIPLE) {
            ++num_triple;
        } else if (bond_orders[i] == chemfiles::Bond::AROMATIC) {
            ++num_aromatic;
        } else if (bond_orders[i] == chemfiles::Bond::AMIDE) { // Strong evidence
            return N_am;
        }
    }

    auto numnonmetal = num_nonmetal(mol_, atom);

    if (numnonmetal >= 4 && num_double == 0 && num_triple == 0) {
        return N_4;
    }

    if (num_aromatic == 2) {
        return N_ar;
    }

    if (numnonmetal == 1 && num_triple > 0) {
        return N_1;
    }

    if (numnonmetal == 2 && (num_double == 2 || num_triple == 1) ) {
        return N_1;
    }

    if (num_amide > 0 && num_deloc != 0) { // Changed! original checks numnonmetal as well
        return N_am;
    }

    if (numnonmetal == 3) {
        if (num_double != 0 || num_triple != 0 || num_aromatic != 0 ||
            num_deloc != 0) {
            return N_pl3;
        }

        return N_3;
    }

    return N_2;
}

size_t Sybyl::assign_carbon_3d_(AtomVertex& atom) {
    auto num_neighbors = atom.neighbor_count();
    if (num_neighbors >= 4) {
        return C_3;
    } else if (num_neighbors == 1) {
        auto dist = mol_.frame().distance(atom, atom[0]);
        if (dist > 1.41)
            return C_3;
        if (dist <= 1.22)
            return C_1;    
        return C_2;
    } else {
        auto avgAngle = 0.0;
        size_t angCount = 0;
        for (size_t n1 = 0; n1 < atom.neighbor_count(); ++n1) {
            for (size_t n2 = n1 + 1; n2 < atom.neighbor_count(); ++n2) {
                avgAngle += mol_.frame().angle(atom[n1], atom, atom[n2]);
                ++angCount;
            }
        }

        avgAngle /= static_cast<double>(angCount);
        avgAngle *= 180.0 / 3.14149;

        if (avgAngle > 160.0) {
            return C_1;
        }

        if (avgAngle <= 115.0) {
            return C_3;
        }

        return C_2;
    }
}

size_t Sybyl::assign_nitrogen_3d_(AtomVertex& atom) {
    auto numnonmetal = num_nonmetal(mol_, atom);
    if (numnonmetal == 4) {
        return N_4;
    } else if (numnonmetal == 1) {
        auto dist = mol_.frame().distance(atom, atom[0]);
        if (dist > 1.2)
            return N_3;
        return N_1;
    } else if (numnonmetal == 3) {
        for (auto neighbor : atom) {
            auto bondee = mol_[neighbor];
            if (bondee.atomic_number() == 6 &&
                freeOxygens(mol_, bondee) == 1) {
                return N_am;
            }
        }

        auto sumAngle = 0.0;
        sumAngle += mol_.frame().angle(atom[0], atom, atom[1]);
        sumAngle += mol_.frame().angle(atom[0], atom, atom[2]);
        sumAngle += mol_.frame().angle(atom[1], atom, atom[2]);
        sumAngle *= 180 / 3.14149;

        if (sumAngle >= 350.0) {
            return N_pl3;
        }

        return N_3;
    }

    auto avgAngle = 0.0;
    size_t angCount = 0;
    for (size_t n1 = 0; n1 < atom.neighbor_count(); ++n1) {
        for (size_t n2 = n1 + 1; n2 < atom.neighbor_count(); ++n2) {
            avgAngle += mol_.frame().angle(atom[n1], atom, atom[n2]);
            ++angCount;
        }
    }

    avgAngle /= static_cast<double>(angCount);
    avgAngle *= 180.0 / 3.14149;

    if (avgAngle > 160.0) {
        return N_1;
    }

    return N_2;
}

size_t Sybyl::assign_oxygen_(AtomVertex& atom) {
    auto numnonmetal = num_nonmetal(mol_, atom);
    if (numnonmetal == 1) {
        auto bondee = atom[0];
        size_t freeOxy = freeOxygens(mol_, atom);

        if (bondee.atomic_number() == 6 &&
            bondee.neighbor_count() == 3 && freeOxy >= 2) {
            return O_co2;
        }

        if (bondee.atomic_number() == 15 && freeOxy >= 2) {
            return O_co2;
        }
    }

    if (numnonmetal == 0 || atom.neighbor_count() >= 2) {
        return O_3;
    }

    return O_2;
}

size_t Sybyl::assign_sulfur_(AtomVertex& atom) {
    auto numnonmetal = num_nonmetal(mol_, atom);

    if (numnonmetal == 3 && freeOxygens(mol_, atom) == 1) {
        return S_o;
    }

    if (numnonmetal == 4 && freeOxygens(mol_, atom) == 2) {
        return S_o2;
    }

    if (atom.neighbor_count() == 2) {
        return S_3;
    } 

    return S_2;
}

bool Sybyl::is_aromatic(size_t atom_id) const {
    switch (atom_types_[atom_id]) {
        case sybyl::C_ar: return true; break;
        case sybyl::N_ar: return true; break;
        default: return false; break;
    } //TODO: Maybe add the ability to check for presence in ring?
}

Hybridization Sybyl::hybridization(size_t atom_id) const {
    return Hybridization::UNKNOWN;
}

template<> std::string Spear::atomtype_name_for_id<Sybyl>(size_t id) {
    return sybyl_unmask[id];
}

template<> size_t Spear::atomtype_id_for_name<Sybyl>(std::string name) {
    return sybyl_mask.at(name);
}

template<> size_t atomtype_id_count<Sybyl>() {
    assert(sybyl_mask.size() == sybyl::Du + 1);
    assert(sybyl_mask.size() == sizeof(sybyl_unmask) / sizeof(char*));
    return sybyl_mask.size();
}
