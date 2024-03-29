// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include "spear/atomtypes/Sybyl.hpp"

#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Geometry.hpp"

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
    Tl, Tm, U,  V,  W,  Xe, Y,  Yb, Zn, Zr, Du
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

static size_t num_nonmetal(const AtomVertex& atom) {
    auto neighs = atom.neighbors();
    auto counts = std::count_if(
        neighs.begin(), neighs.end(),
        [](AtomVertex n){ return n.is_non_metal();}
    );
    return static_cast<std::size_t>(counts);
}

static size_t freeOxygens(const Spear::AtomVertex& atom) {
    size_t freeOxygens = 0;
    for (auto bondee : atom.neighbors()) {
        if (bondee.atomic_number() == 8 && num_nonmetal(bondee) == 1) {
            ++freeOxygens;
        }
    }

    return freeOxygens;
}

static bool is_decloc(const Spear::AtomVertex& atom) {
    for (auto bond : atom.bonds()) {
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

Sybyl::Sybyl(const Molecule& mol, TypingMode mode) :
    mol_(mol) {

    super::resize(mol.size(), sybyl::Du);
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
        case Element::C:
            set(atom, assign_carbon_3d_(atom));
            break;
        case Element::N:
            set(atom, assign_nitrogen_3d_(atom));
            break;
        default:
            break;
        }
    }
}

void Sybyl::type_atoms_topo_(const AtomVertex& atom) {
    switch (atom.atomic_number()) {
    case Element::H:
        set(atom, sybyl::H);
        break;
    case Element::C:
        set(atom, assign_carbon_topo_(atom));
        break;
    case Element::N:
        set(atom, assign_nitrogen_topo_(atom));
        break;
    case Element::O:
        set(atom, assign_oxygen_(atom));
        break;
    case Element::S:
        set(atom, assign_sulfur_(atom));
        break;
    case Element::P:
        set(atom, sybyl::P_3);
        break;
    case Element::Co:
        set(atom, sybyl::Co_oh);
        break;
    case Element::Ru:
        set(atom, sybyl::Ru_oh);
        break;
    case Element::Ti:
        set(atom, atom.degree() <= 4 ? Ti_th
                                     : Ti_oh);
        break;
    case Element::Cr:
        set(atom, atom.degree() <= 4 ? Cr_th
                                     : Cr_oh);
        break;
    default:
        set(atom, sybyl_mask.at(Element::Name[atom.atomic_number()]));
        break;
    }
}

void Sybyl::type_atoms_topo_() {
    for (auto atom : mol_) {
        type_atoms_topo_(atom);
    }
}

size_t Sybyl::assign_carbon_topo_(const AtomVertex& atom){
    size_t num_double = 0, num_triple = 0, num_aromatic = 0;
    size_t num_nitrogen = 0;

    for (auto bond : atom.bonds()) {
        if (bond.source().atomic_number() == Element::N) {
            num_nitrogen += static_cast<size_t>(freeOxygens(bond.source()) == 0);
        }

        if (bond.target().atomic_number() == Element::N) {
            num_nitrogen += static_cast<size_t>(freeOxygens(bond.target()) == 0);
        }

        if (bond.order() == Bond::DOUBLE) {
            ++num_double;
        } else if (bond.order() == Bond::TRIPLE) {
            ++num_triple;
        } else if (bond.order() == Bond::AROMATIC) {
            ++num_aromatic;
        }
    }

    auto degree = atom.degree();

    if (degree >= 4 && num_double == 0 && num_triple == 0) {
        return C_3;
    }

    if (num_aromatic >= 2) {
        return C_ar;
    }

    if (num_nitrogen == 3 && degree == 3) {
        // TODO Need to check for acyclic...
        return C_cat;
    }

    if (num_triple == 1 || num_double == 2) {
        return C_1;
    }

    return C_2;
}

size_t Sybyl::assign_nitrogen_topo_(const AtomVertex& atom) {
    size_t num_double = 0, num_triple = 0, num_aromatic = 0;
    size_t num_amide = 0, num_deloc = 0, num_nitrogens = 0;

    auto update = [&] (const AtomVertex& be) {
        if (be.atomic_number() == 6) {
            num_amide += static_cast<size_t>(freeOxygens(be) != 0);
            num_deloc += static_cast<size_t>(is_decloc(be));
            num_nitrogens = 0;
            for (auto bondee_bondee : be.neighbors()) {
                num_nitrogens += static_cast<size_t>(bondee_bondee.atomic_number() == Element::N);
            }
        }
    };

    for (auto bond : atom.bonds()) {
        update(bond.source());
        if (num_nitrogens >= 3 && bond.order() == Bond::DOUBLE) {
            return N_pl3;
        }  

        update(bond.target());
        if (num_nitrogens >= 3 && bond.order() == Bond::DOUBLE) {
            return N_pl3;
        }

        if (bond.order() == Bond::DOUBLE) {
            ++num_double;
        } else if (bond.order() == Bond::TRIPLE) {
            ++num_triple;
        } else if (bond.order() == Bond::AROMATIC) {
            ++num_aromatic;
        } else if (bond.order() == Bond::AMIDE) { // Strong evidence
            return N_am;
        }
    }

    auto numnonmetal = num_nonmetal(atom);

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

size_t Sybyl::assign_carbon_3d_(const AtomVertex& atom) {
    auto num_neighbors = atom.degree();
    if (num_neighbors >= 4) {
        return C_3;
    } else if (num_neighbors == 1) {
        auto dist = distance(atom.position(), atom[0].position());
        if (dist > 1.41)
            return C_3;
        if (dist <= 1.22)
            return C_1;    
        return C_2;
    } else {
        auto avgAngle = 0.0;
        size_t angCount = 0;
        for (size_t n1 = 0; n1 < atom.degree(); ++n1) {
            for (size_t n2 = n1 + 1; n2 < atom.degree(); ++n2) {
                avgAngle += angle(atom[n1].position(),
                                  atom.position(),
                                  atom[n2].position());
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

size_t Sybyl::assign_nitrogen_3d_(const AtomVertex& atom) {
    auto numnonmetal = num_nonmetal(atom);
    if (numnonmetal == 4) {
        return N_4;
    } else if (numnonmetal == 1) {
        auto dist = distance(atom.position(), atom[0].position());
        if (dist > 1.2)
            return N_3;
        return N_1;
    } else if (numnonmetal == 3) {
        for (auto bondee : atom.neighbors()) {
            if (bondee.atomic_number() == 6 &&
                freeOxygens(bondee) == 1) {
                return N_am;
            }
        }

        auto sumAngle = 0.0;
        sumAngle += angle(atom[0].position(), atom.position(), atom[1].position());
        sumAngle += angle(atom[0].position(), atom.position(), atom[2].position());
        sumAngle += angle(atom[1].position(), atom.position(), atom[2].position());
        sumAngle *= 180 / 3.14149;

        if (sumAngle >= 350.0) {
            return N_pl3;
        }

        return N_3;
    }

    auto avgAngle = 0.0;
    size_t angCount = 0;
    for (size_t n1 = 0; n1 < atom.degree(); ++n1) {
        for (size_t n2 = n1 + 1; n2 < atom.degree(); ++n2) {
            avgAngle += angle(atom[n1].position(), atom.position(), atom[n2].position());
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

size_t Sybyl::assign_oxygen_(const AtomVertex& atom) {
    auto numnonmetal = num_nonmetal(atom);
    if (numnonmetal == 1) {
        auto bondee = atom[0];
        size_t freeOxy = freeOxygens(atom);

        if (bondee.atomic_number() == 6 &&
            bondee.degree() == 3 && freeOxy >= 2) {
            return O_co2;
        }

        if (bondee.atomic_number() == 15 && freeOxy >= 2) {
            return O_co2;
        }
    }

    if (numnonmetal == 0 || atom.degree() >= 2) {
        return O_3;
    }

    return O_2;
}

size_t Sybyl::assign_sulfur_(const AtomVertex& atom) {
    auto numnonmetal = num_nonmetal(atom);

    if (numnonmetal == 3 && freeOxygens( atom) == 1) {
        return S_o;
    }

    if (numnonmetal == 4 && freeOxygens(atom) == 2) {
        return S_o2;
    }

    if (atom.degree() == 2) {
        return S_3;
    } 

    return S_2;
}

bool Sybyl::is_aromatic(size_t atom_id) const {
    switch (get(atom_id)) {
    case sybyl::C_ar:
    case sybyl::N_ar:
        return true;
    default: return false;
    } //TODO: Maybe add the ability to check for presence in ring? Oxygen aromaticity?
}

bool Sybyl::is_planar(size_t atom_id) const {
    switch (get(atom_id)) {
    // Carbon
    case sybyl::C_ar: case sybyl::C_2: case sybyl::C_1: case sybyl::C_cat:
    // Nitrogen
    case sybyl::N_ar: case sybyl::N_2: case sybyl::N_1: case sybyl::N_am:
    case sybyl::N_pl3:
    // Oxygen
    case sybyl::O_co2: case sybyl::O_2:
    // Sulfur
    case sybyl::S_2:
        return true;
    // Not planar
    default:
        return false;
    }
}

Hybridization Sybyl::hybridization(size_t atom_id) const {
    switch (get(atom_id)) {
    case sybyl::C_1: case sybyl::N_1:
        return Hybridization::SP;

    case sybyl::C_2:  case sybyl::N_2:  case sybyl::O_2:   case sybyl::S_2:
    case sybyl::C_ar: case sybyl::N_ar: case sybyl::O_co2:
    case sybyl::C_cat:
        return Hybridization::SP2;

    case sybyl::C_3:   case sybyl::N_3: case sybyl::N_4:   case sybyl::O_3:
    case sybyl::S_3:   case sybyl::P_3: case sybyl::Ti_th: case sybyl::Cr_th:
    case sybyl::N_pl3: case sybyl::S_o2: case sybyl::S_o:  case sybyl::N_am:
        return Hybridization::SP3;

    default:
        return Hybridization::FORCED;
    }
}

size_t Sybyl::add_atom(size_t idx) {
    switch(mol_[idx].atomic_number()) {
    case Element::H:
        push_back(sybyl::H);
        break;
    case Element::C:
        push_back(sybyl::C_3);
        break;
    case Element::N:
        push_back(sybyl::N_3);
        break;
    case Element::O:
        push_back(sybyl::O_3);
        break;
    case Element::S:
        push_back(sybyl::S_3);
        break;
    case Element::P:
        push_back(sybyl::P_3);
        break;
    default:
        push_back(sybyl_mask.at(Element::Name[mol_[idx].atomic_number()]));
        break;
    }

    return back();

}

void Sybyl::add_bond(size_t idx1, size_t idx2, Bond::Order /*unused*/) {
    type_atoms_topo_(mol_[idx1]);
    type_atoms_topo_(mol_[idx2]);
}

template<> std::string Spear::atomtype_name_for_id<Sybyl>(size_t id) {
    return sybyl_unmask[id];
}

template<> size_t Spear::atomtype_id_for_name<Sybyl>(const std::string& name) {
    return sybyl_mask.at(name);
}

template<> size_t Spear::atomtype_id_count<Sybyl>() {
    assert(sybyl_mask.size() == sybyl::Du + 1);
    assert(sybyl_mask.size() == sizeof(sybyl_unmask) / sizeof(char*));
    return sybyl_mask.size();
}
