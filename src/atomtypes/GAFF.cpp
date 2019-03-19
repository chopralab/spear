#include "spear/atomtypes/GAFF.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/FunctionalGroup.hpp"

#include <map>
#include <string>
#include <locale>

using namespace Spear;

enum gaff : uint64_t {
    X = 0, DU = 1,
    hc, ha, ho, hn, hs, hp, h1, h2, h3, h4, h5,
    c,  ca, c1, c2, c3, cg, cz, ce, cc, cu, cv, cd, cp, cx, cy,
    f,  cl, br, i,
    p2, p3, pc, pe, pb, px, p4, p5, py,
    n,  n1, n2, n3, n4, na, nb, nc, nd, ne, nh, no,
    o,  os, oh, ow, 
    s,  ss, s2, sh, s4, sx, sy, s6, 
    He, Li, Be, B,  Ne, Na, Mg, Al, Si, Ar,
    K,  Ca, Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Kr,
    Rb, Sr, Y,  Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, Xe,
    Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
        Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
};

GAFF::GAFF(const Molecule& mol) : mol_(mol) {
    atom_types_.reserve(mol.size());

    for (auto av : mol) {
        add_atom(av);
    }
}

bool GAFF::is_aromatic(size_t atom_id) const {
}

bool GAFF::is_planar(size_t atom_id) const {
}

static bool follows_rule(const AtomVertex& vert, const std::string& rule) {
    FunctionalGroup fg(rule);
    auto groups = find_functional_groups(vert.br(), fg);

    for (auto match : groups) {
        if (match.size() && match[0] == vert) {
            return true;
        }
    }

    return false;
}

static size_t type_oxygen(const AtomVertex& av) {
    auto sssrs = av.sssrs();
    auto sssrs_count = std::distance(sssrs.first, sssrs.second);
    switch(av.degree() + av.implicit_hydrogens()) {
    case 1: // Must be SP2
        return gaff::o;
    case 2:
        // Between two carbonyls (acid anhydrous)
        if (follows_rule(av, "[OX2]([#6X3]=[OX1])[#6X3]=[OX1]")) return gaff::os;
        // Aliphatic ester
        if (follows_rule(av, "[OX2]([#6X3]=[OX1])[#6X4]"))       return gaff::os;

        if (av.total_hydrogens() == 0 && sssrs_count == 0) {
            if (follows_rule(av, "[OX2]([#6X3])([CX3]=[OX1])"))  return gaff::os;
            if (follows_rule(av, "[OX2]([#6X4])([NX3]=[OX1])"))  return gaff::os;
            if (follows_rule(av, "[OX2][NX3A]=[OX1]"))           return gaff::o;
            if (follows_rule(av, "[OX2][NX2A]=[OX1]"))           return gaff::o;
            if (follows_rule(av, "[OX2][CX3A]=[OX1]"))           return gaff::o;
        }

        if (av.total_hydrogens() == 2) return gaff::ow;
        if (av.total_hydrogens() == 1) return gaff::oh;

        if (follows_rule(av, "[OX2&H0]([#16X4])[#6]"))    return gaff::os;
        if (follows_rule(av, "[OX2r5]([#15X3][Cl])[#6]")) return gaff::os;
        if (follows_rule(av, "[OX2&H0]([#15X4])[#6]"))    return gaff::os;

        if (follows_rule(av, "[OX2&H0][#16X4]-[#6X4]"))   return gaff::o;
        if (follows_rule(av, "[OX2&H0][#16X4]=[OX1]"))    return gaff::o;
        if (follows_rule(av, "[OX2&H0][#15X4]=[OX1]"))    return gaff::o;

        return gaff::os; // Give up, guess ester / ether
    case 3: // This is odd, but check if there's a hydrogen just in case
        return av.total_hydrogens() >= 1? gaff::oh : gaff::os;
    default: // No idea what to do
        return gaff::os;
    }
}

static size_t type_phosphorus(const AtomVertex& av) {
    auto sssrs = av.sssrs();
    auto sssrs_count = std::distance(sssrs.first, sssrs.second);

    switch(av.degree() + av.implicit_hydrogens()) {
    case 1:
        return gaff::p2;
    case 2:
        if (sssrs_count && av.is_aromatic()) return gaff::pb;
        if (sssrs_count) {
            if (follows_rule(av, "[#15X2](#*)-*#*"))    return gaff::pc;
            if (follows_rule(av, "[#15X2](=[#6])[#8]")) return gaff::pc;
            if (follows_rule(av, "[#15X2](=*)-*#*"))    return gaff::pc;
            if (follows_rule(av, "[#15X2](=*)-*=*"))    return gaff::pc;
        } else {
            if (follows_rule(av, "[#15X2](#*)-*#*"))    return gaff::pe;
            if (follows_rule(av, "[#15X2](=*)-*#*"))    return gaff::pe;
            if (follows_rule(av, "[#15X2](=*)-*=*"))    return gaff::pe;
        }
        return gaff::p2;
    case 3:
        if (follows_rule(av, "[#15X3]=O"))                 return gaff::p4;
        if (follows_rule(av, "[#15X3]=S"))                 return gaff::p4;
        if (follows_rule(av, "[#15X3](=[#8])(=[#8])[#8]")) return gaff::p4;
        if (follows_rule(av, "[#15X3](=*)-*#*"))           return gaff::px;
        if (follows_rule(av, "[#15X3](=*)-*=*"))           return gaff::px;
        return gaff::p3;
    case 4:
        if (follows_rule(av, "[#15X4](=*)-*#*")) return gaff::py;
        if (follows_rule(av, "[#15X4](=*)-*=*")) return gaff::py;
        return gaff::p5;
    case 5:
    case 6:
        return gaff::p5;
    default: // no idea
        return gaff::p5;
    }
}

static size_t type_sulfur(const AtomVertex& av) {
    switch(av.degree() + av.implicit_hydrogens()) {
    case 1:
        return gaff::s;
    case 2:
        if (av.total_hydrogens())          return gaff::sh;
        if (follows_rule(av, "[SX2]=*"))   return gaff::s2;
        if (follows_rule(av, "[SX2]#*"))   return gaff::s2;
        return gaff::ss;
    case 3:
        if (follows_rule(av, "[SX3](=*)-*=*"))   return gaff::sx;
        if (follows_rule(av, "[SX3](=*)-*#*"))   return gaff::sx;
        return gaff::s4;
    case 4: // check for sulfoxides and sulfates
        if (follows_rule(av, "[SX4](=*)-*=*"))   return gaff::sy;
        if (follows_rule(av, "[SX4](=*)-*#*"))   return gaff::sy;
        // fall through
    case 5:
    case 6:
        return gaff::s6;
    }
}

size_t GAFF::add_atom(size_t new_idx) {
    auto av = mol_[new_idx];

    // NEED TO RESIZE FIRST

    switch(av.atomic_number()) {
    case Element::O:
        atom_types_[new_idx] = type_oxygen(mol_[new_idx]);
        break;
    case Element::P:
        atom_types_[new_idx] = type_phosphorus(mol_[new_idx]);
        break;
    case Element::S:
        atom_types_[new_idx] = type_sulfur(mol_[new_idx]);
        break;
    case Element::F:
        atom_types_[new_idx] = gaff::f;
        break;
    case Element::Cl:
        atom_types_[new_idx] = gaff::cl;
        break;
    case Element::Br:
        atom_types_[new_idx] = gaff::br;
        break;
    case Element::I:
        atom_types_[new_idx] = gaff::i;
        break;
    }

    return atom_types_[new_idx];
}

void GAFF::add_bond(size_t idx1, size_t idx2, Bond::Order bo) {
}

template<> std::string atomtype_name_for_id<GAFF>(size_t id) {
}

template<> size_t atomtype_id_for_name<GAFF>(std::string name) {
}

template<> size_t atomtype_id_count<GAFF>() {
}

template<> double van_der_waals<GAFF>(size_t id) {
}
