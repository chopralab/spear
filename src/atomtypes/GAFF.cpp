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
    hc, ha, ho, hn, hs, hp, h1, h2, h3, h4, h5, hw,
        hx, // Unused
    c,  ca, c1, c2, c3, cg, cz, ce, cc, cu, cv, cd, cp, cx, cy,
        cs, cq, cf, // Unused
    f,  cl, br, i,
    p2, p3, pc, pe, pb, px, p4, p5, py,
        pd, pf, // Unused
    n,  n1, n2, n3, n4, na, nb, nc, nd, ne, nh, no,
        nf, ns, nt, nx, ny, nz, nP, nu, nv, n7, n8, n9, // Unused
    o,  os, oh, ow, 
    s,  ss, s2, sh, s4, sx, sy, s6, 
    He, Li, Be, B,  Ne, Na, Mg, Al, Si, Ar,
    K,  Ca, Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Kr,
    Rb, Sr, Y,  Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, Xe,
    Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
        Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
};

const char* const gaff_unmask[]{
    "X", "DU",
    "hc", "ha", "ho", "hn", "hs", "hp", "h1", "h2", "h3", "h4", "h5", "hw", "hx",
    "c",  "ca", "c1", "c2", "c3", "cg", "cz", "ce", "cc", "cu", "cv", "cd", "cp", "cx", "cy", "cs", "cq", "cf",
    "f",  "cl", "br", "i",
    "p2", "p3", "pc", "pe", "pb", "px", "p4", "p5", "py", "pd", "pf",
    "n",  "n1", "n2", "n3", "n4", "na", "nb", "nc", "nd", "ne", "nh", "no",
        "nf", "ns", "nt", "nx", "ny", "nz", "n+", "nu", "nv", "n7", "n8", "n9",
    "o",  "os", "oh", "ow", 
    "s",  "ss", "s2", "sh", "s4", "sx", "sy", "s6", 
    "He", "Li", "Be", "B",  "Ne", "Na", "Mg", "Al", "Si", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
          "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
};

#define AMINO_ACID_CORE {"OXT","o"},{"CA","c3"},{"C","c"},{"O","o"},{"N","n"},{"HA","h4"}

const std::map<std::string, std::map<std::string,std::string>> standard_residues {
    {"ALA",{AMINO_ACID_CORE,{"CB","c3"}}},
    {"ARG",{AMINO_ACID_CORE,{"CB","c3"},{"CD","c3"},{"CG","c3"},{"CZ","cz"},{"NE","nh"},{"NH1","nh"},{"NH2","nh"}}},
    {"ASN",{AMINO_ACID_CORE,{"CB","c3"},{"CG","c"},{"ND2","n"},{"OD1","o"}}},
    {"ASP",{AMINO_ACID_CORE,{"CB","c3"},{"CG","c"},{"OD1","o"},{"OD2","o"}}},
    {"CYS",{AMINO_ACID_CORE,{"CB","c3"},{"SG","s"}}},
    {"GLN",{AMINO_ACID_CORE,{"CB","c3"},{"CD","c"},{"CG","c3"},{"NE2","n"},{"OE1","o"}}},
    {"GLU",{AMINO_ACID_CORE,{"CB","c3"},{"CD","c"},{"CG","c3"},{"OE1","o"},{"OE2","o"}}},
    {"GLY",{AMINO_ACID_CORE}},
    {"HIS",{AMINO_ACID_CORE,{"CB","c3"},{"CD2","cc"},{"CE1","cc"},{"CG","cc"},{"ND1","nc"},{"NE2","nc"}}},
    {"ILE",{AMINO_ACID_CORE,{"CB","c3"},{"CD1","c3"},{"CG1","c3"},{"CG2","c3"},}},
    {"LEU",{AMINO_ACID_CORE,{"CB","c3"},{"CD1","c3"},{"CD2","c3"},{"CG","c3"},}},
    {"LYS",{AMINO_ACID_CORE,{"CB","c3"},{"CD","c3"},{"CE","c3"},{"CG","c3"},{"NZ","n3"},}},
    {"MET",{AMINO_ACID_CORE,{"CB","c3"},{"CE","c3"},{"CG","c3"},{"SD","s"}}},
    {"MSE",{AMINO_ACID_CORE,{"CB","c3"},{"CE","c3"},{"CG","c3"},{"SE","DU"}}},
    {"PHE",{AMINO_ACID_CORE,{"CB","c3"},{"CD1","ca"},{"CD2","ca"},{"CE1","ca"},{"CE2","ca"},{"CG","ca"},{"CZ","ca"}}},
    {"PRO",{AMINO_ACID_CORE,{"CB","c3"},{"CD","c3"},{"CG","c3"}}},
    {"SER",{AMINO_ACID_CORE,{"CB","c3"},{"OG","oh"}}},
    {"THR",{AMINO_ACID_CORE,{"CB","c3"},{"CG2","c3"},{"OG1","oh"}}},
    {"TRP",{AMINO_ACID_CORE,{"CB","c3"},{"CD1","ca"},{"CD2","ca"},{"CE2","ca"},{"CE3","ca"},{"CG","ca"},{"CH2","ca"},{"CZ2","ca"},{"CZ3","ca"},{"NE1","na"}}},
    {"TYR",{AMINO_ACID_CORE,{"CB","c3"},{"CD1","ca"},{"CD2","ca"},{"CE1","ca"},{"CE2","ca"},{"CG","ca"},{"CZ","ca"},{"OH","oh"}}},
    {"VAL",{AMINO_ACID_CORE,{"CB","c3"},{"CG1","c3"},{"CG2","c3"}}},
    {"HOH",{{"O","ow"}}},
};

GAFF::GAFF(const Molecule& mol) : mol_(mol) {
    reserve(mol.size());

    for (auto av : mol) {
        add_atom_(av);
    }
}

bool GAFF::is_aromatic(size_t /*atom_id*/) const {
	return false;
}

bool GAFF::is_planar(size_t /*atom_id*/) const {
	return false;
}

static bool follows_rule(const AtomVertex& vert, const std::string& rule) {
    FunctionalGroup fg(rule);
    auto groups = find_functional_groups(vert.br(), fg);

    for (auto match : groups) {
        if (match.size() != 0 && match[0] == vert) {
            return true;
        }
    }

    return false;
}

static size_t num_withdrawing(const AtomVertex& av) {
    size_t count = 0;
    for (auto neighbor : av.neighbors()) {
        switch (neighbor.atomic_number()) {
        case Element::N:
        case Element::O:
        case Element::F:
        case Element::S:
        case Element::Cl:
        case Element::Br:
        case Element::I:
            ++count;
            break;
        default:
            break;
        }
    }
    return count;
}

static size_t type_hydrogen(const AtomVertex& av) {
    if (av.degree() != 1) {
        return gaff::hc; // hc is 'other' carbon
    }

    auto neighbor = av[0]; // The only neighbor
    switch (neighbor.atomic_number()) {
    case Element::O:
        if (neighbor.explicit_hydrogens() >= 2)
            return gaff::hw;
        return gaff::ho;
    case Element::N:
        return gaff::hn;
    case Element::S:
        return gaff::hs;
    case Element::P:
        return gaff::hp;
    case Element::Ge:
        return gaff::DU;
    case Element::C:
        break;
    default:
        return gaff::ha; // ha is also 'other' carbon
    }

    // We know it must be next to a carbon!
    if (neighbor.is_aromatic()) {
        if (num_withdrawing(neighbor) == 1)        return gaff::h4;
        if (num_withdrawing(neighbor) == 2)        return gaff::h5;
        if (follows_rule(av, "[#1X1][c]-[#15X3]")) return gaff::h4;
        if (follows_rule(av, "[#1X1][c]:[#15X3]")) return gaff::h4;
        return gaff::ha;
    }

    switch (num_withdrawing(neighbor)) {
    case 1:
        return gaff::h1;
    case 2:
        return gaff::h2;
    case 3:
        return gaff::h3;
    default:
        break;
    }

    if (neighbor.expected_bonds() == 4) {
        return gaff::hc;
    }

    return gaff::ha; // Must be conjugated somehow
}

static size_t type_carbon(const AtomVertex& av) {
    auto sssrs = av.sssrs();
    auto sssrs_count = std::distance(sssrs.first, sssrs.second);
    std::set<size_t> r;
    for (auto riter = sssrs.first; riter != sssrs.second; ++riter) {
        r.insert(riter->second.size()); // ring size
    }
    switch(av.degree() + av.implicit_hydrogens()) {
    case 1:
        return gaff::c1;
    case 2:
        if (follows_rule(av, "[#6](#[#6])[#1]"))    return gaff::c1;
        if (follows_rule(av, "[#6](#[#7])[#6X3]"))  return gaff::cg;
        if (follows_rule(av, "[#6](#*)-*#*"))       return gaff::cg;
        if (follows_rule(av, "[#6](#*)-*=*"))       return gaff::cg;
        if (follows_rule(av, "[#6](=*)-*=*"))       return gaff::cg;
        return gaff::c1;
    case 3:
        if(r.count(8) != 0) {
            if (follows_rule(av, "[#6]=*-*=*"))       return gaff::cc;
            if (follows_rule(av, "[#6](=[#6])([#6])[#1]"))      return gaff::c2;
            if (follows_rule(av, "[#6](=[#6]-*)(-[#6]=*)[#1]")) return gaff::cc;
        }
        if(r.count(7) != 0) {
            if (follows_rule(av, "[#6]([#6X3]:[#6X3][#6X3]=[O])[#6X3][#1]")) return gaff::cc;
            if (follows_rule(av, "[#6](:[#6X3][#6X3]=[O])[#6X3][#1]"))       return gaff::cc;
            if (follows_rule(av, "[#6](:[#6X3])([#6X3]=[O])[#1]"))           return gaff::cc;
            if (follows_rule(av, "[#6](=[#6][#15])([#6])[#1]"))              return gaff::cc;
            if (follows_rule(av, "[#6](!:[#6])([#6]=[#6])"))                 return gaff::cc;
            if (follows_rule(av, "[#6](=[#6X3])([#15X3])[#1]"))              return gaff::cc;
            if (follows_rule(av, "[#6](=[#6][#15])([#6])[#1]"))              return gaff::cc;
            if (follows_rule(av, "[#6][#6]=[#6]-[#15]"))                     return gaff::cc;
            if (follows_rule(av, "[#6](=[#6])([#6])[#1]"))                   return gaff::c2;
        }
        if(r.count(6) != 0) {
            if (follows_rule(av, "c([nX3])([cX3][cX3]=O)"))                       return gaff::cc;
            if (follows_rule(av, "c([cX3]=O)"))                                   return gaff::cd;
            if (follows_rule(av, "c(:[c][cX3][NX3A])([nX3])"))                    return gaff::cc;
            if (follows_rule(av, "c(:[c][nX3])([c][NX3A])"))                      return gaff::cd;
            if (follows_rule(av, "c(:nc=[OA])([NX3A])c"))                         return gaff::cd;
            if (follows_rule(av, "c=O"))                                          return gaff::c;
            if (follows_rule(av, "ccn"))                                          return gaff::ca;
            if (follows_rule(av, "cn"))                                           return gaff::ca;
            if (follows_rule(av, "[#6](=[#6X3][#6X3](=[O])[O])([#6X3]=[O])[#1]")) return gaff::cc;
            if (follows_rule(av, "[#6](=[#6X3])([#6X3](=[O])[O])[#6X3]=[O]"))     return gaff::cc;
            if (follows_rule(av, "[#6](=[#6X3][#7])([#6X3]=[#7])[#1]"))           return gaff::cc;
            if (follows_rule(av, "[#6]([#6X4]-[#6X3])(=[#6X3]-[#8])[#1]"))        return gaff::c2;
            if (follows_rule(av, "[#6](-[#6X4]-[#6X4])(=[#6X3]-[#8])[#1]"))       return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3]-[F])([#6X4]-[#6X4])[#1]"))         return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3][#6X4])([#6X3])[#1]"))              return gaff::c2;
            if (follows_rule(av, "[#6]([#6X4])(=[#6X3]-[#6X3])[#1]"))             return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3]-[#6X4])([#6X3]=[#8])[#1]"))        return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3]-[#6X4])([#8]-[#6X3])[#1]"))        return gaff::c2;
            if (follows_rule(av, "[#6](=[#6]-[#6X4])([#8]-[#6X4])[#1]"))          return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3]-[#6X3]=[#8X1])([#6X3])[#1]"))      return gaff::cd;
            //if (follows_rule(av, "[#6](:[#6])([#6]@[#7])[#1]"))                   return gaff::cc; Needs to be fixed
            if (follows_rule(av, "[#6](=[#6X3])([#6X3])[F]"))                     return gaff::ca;
            if (follows_rule(av, "[#6](=[#6X3])([#6X3])[Cl]"))                    return gaff::ca;
            if (follows_rule(av, "[#6](=[#6X3])([#6X3])[Br]"))                    return gaff::ca;
            if (follows_rule(av, "[#6](=[#6X3])([#6X3])[I]"))                     return gaff::ca;
            if (follows_rule(av, "[#6](=[#6X3])([#6X4])[#6X4]"))                  return gaff::c2;
            if (follows_rule(av, "[#6](=[#6])([#6X4])[#1]"))                      return gaff::c2;
        }

        if (follows_rule(av, "c-c"))                                              return gaff::cp;

        if (r.count(5) != 0) {
            if (follows_rule(av, "[#6]=O"))                                       return gaff::c;
            if (follows_rule(av, "[#6]=S"))                                       return gaff::c;
            if (follows_rule(av, "[#6](=[#6X3]-[#16X2])(-[#6X4]-[#6X4])[#1]"))    return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3]-[#6X4])([#6X4]-[#6X4])[#1]"))      return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3][#8X2])([#6X4]-[#6X3])[#1]"))       return gaff::c2;
            if (follows_rule(av, "[#6](-[#6X3]=[#6X3])(=[#6X3]-[#6X4])[#1]"))     return gaff::cd;
            if (follows_rule(av, "[#6]([#6X4])(=[#6X3]-[#8])[#1]"))               return gaff::c2;
            if (follows_rule(av, "[#6](=[#6X3])(-[#6X4]-[#8])[#1]"))              return gaff::c2;
            if (follows_rule(av, "[#6](-[#6X4])(=[#6X3]-[#6X3])[#1]"))            return gaff::cc;
            if (follows_rule(av, "[#6](=[#6X3][#6X3]=[#8X1])[#6X3]"))             return gaff::cc;
            if (follows_rule(av, "[#6](=[#6X3][#6X3])[#6X3]=[#8X1]"))             return gaff::cc;
            if (follows_rule(av, "[#6](:[#6X3]-[#16X2])([#6X3])[#1]"))            return gaff::cd;
            if (follows_rule(av, "[#6](:[#6X3])(-[#16X2])[#1]"))                  return gaff::cc;
            if (follows_rule(av, "[#6](=[#6X3])(-[#6X3]=[#8X1])[#1]"))            return gaff::cc;
            if (follows_rule(av, "[#6](=[#6])([#6])[#1]"))                        return gaff::c2;
            if (follows_rule(av, "[#6](-[#6]=[#6]-[#6])"))                        return gaff::cc;
            if (follows_rule(av, "[c][nX3]"))                                     return gaff::cc;
            if (follows_rule(av, "[c][nX2]"))                                     return gaff::cd;
            if (follows_rule(av, "[c]n"))                                         return gaff::cc;
            if (follows_rule(av, "[c]o"))                                         return gaff::cc;
        }
        if (r.count(4) != 0) {
            if (follows_rule(av, "[#6]=O"))                                       return gaff::c;
            if (follows_rule(av, "[#6]-*=*-*"))                                   return gaff::cc;
            if (follows_rule(av, "[#6]=*-*=*"))                                   return gaff::cc;
            return gaff::c;
        }
        if (r.count(3) != 0) {
            if (follows_rule(av, "[#6]=O"))                                       return gaff::c;
            return gaff::cu;
        }

        if (follows_rule(av, "[#6]=[OA]"))                                        return gaff::c;
        if (follows_rule(av, "[#6]=[SA]"))                                        return gaff::c;

        if (sssrs_count == 0) {
            if (follows_rule(av, "[#6H](=*)-*~a"))                                return gaff::ce;
            if (follows_rule(av, "[#6](=*)-*#*"))                                 return gaff::ce;
            if (follows_rule(av, "[#6](=*)-*=*"))                                 return gaff::ce;
            if (follows_rule(av, "[#6](=[#7])([#7])[#7]"))                        return gaff::cz;
        } else {
            if (follows_rule(av, "[#6](=*)-*#*"))                                 return gaff::cc;
            if (follows_rule(av, "[#6](=*)-*=*"))                                 return gaff::ca;
        }

        return av.is_aromatic()? gaff::ca : gaff::c2;
    case 4:
        if (follows_rule(av, "[#6r6]([#6X3])([#6X3])[#6X4]")) return gaff::cc;
        if (follows_rule(av, "[#6r4]"))                       return gaff::cy;
        if (follows_rule(av, "[#6r3]"))                       return gaff::cx;
        return gaff::c3;
    default: // No idea
        return gaff::c3;
    }
}

static size_t type_nitrogen(const AtomVertex& av) {
    auto sssrs = av.sssrs();
    auto sssrs_count = std::distance(sssrs.first, sssrs.second);
    switch(av.degree() + av.implicit_hydrogens()) {
    case 1:
        return gaff::n1;
    case 2:
        if (sssrs_count == 0) {
            if (follows_rule(av, "[#7](#*)-*#*")) return gaff::ne;
            if (follows_rule(av, "[#7](=*)-*#*")) return gaff::ne;
            if (follows_rule(av, "[#7](=*)-*=*")) return gaff::ne;
        } else {
            if (follows_rule(av, "[#7](-[#7X3]-[#6X3])=[#6X3]")) return gaff::nc;
            if (follows_rule(av, "[#7](=*)-*=*"))                return gaff::nc;
            if (follows_rule(av, "[#7](=*)-*#*"))                return gaff::nc;
            if (follows_rule(av, "[#7](#*)-*#*"))                return gaff::nc;
            if (follows_rule(av, "[#7](:[#7])[#7]"))             return gaff::nc;
            if (follows_rule(av, "[#7](:[#7])[#16X2]"))          return gaff::nc;
            if (follows_rule(av, "[#7]([#7])[#7]"))              return gaff::nc;
            if (follows_rule(av, "[#7]=[#6]([#7])[#7]"))         return gaff::nc;

            // special ring cases
            if (follows_rule(av, "[nr5]([#6])[#8,#16]")) return gaff::nc;
            if (follows_rule(av, "[nr6](:[#6X3][NX3])([#6X3]=[O])")) return gaff::nd;
            if (follows_rule(av, "[#7r5](:[#6])[#6]")) return gaff::nb;            

            // SP2 and in a conjugated ring
            if (follows_rule(av, "[#7](=[#7])[#16]")) return gaff::nd;
            if (follows_rule(av, "[#7](:[#6])[#7]"))  return gaff::nd;
            if (follows_rule(av, "[#7]=[#6][#7]"))    return gaff::nd;
            if (follows_rule(av, "[#7]=[#6][#8]"))    return gaff::nd;
            if (follows_rule(av, "[#7]=[#6][#16]"))   return gaff::nd;

            // Label as aromatic (even if they are not technically)
            if (follows_rule(av, "[#7](=*)=*"))    return gaff::nb;
            if (follows_rule(av, "[#7](=*)(-*)"))  return gaff::nb;
        }

        if (av.is_aromatic())                             return gaff::nb;
        if (follows_rule(av, "[#7](=[#7X2])-*"))        return gaff::n2;
        if (follows_rule(av, "[#7](#*)-*"))             return gaff::n1;
        if (follows_rule(av, "[#7](=[#7X2])(=[#7X2])")) return gaff::n1;
        if (follows_rule(av, "[#7](=[#7X2]=[#7X2])"))   return gaff::n1;

        return gaff::n2;
    case 3:
        // Nitro group
        if (follows_rule(av, "[#7H0](=O)=O"))          return gaff::no;
        if (follows_rule(av, "[#7H0](=O)-O"))          return gaff::no;

        // Amides
        if (follows_rule(av, "[#7r6]([#6X3]=[#8])([#6X3])[#6X4]"))      return gaff::n;
        if (follows_rule(av, "[#7r6]([#6X3]=[#8])([#6X3])[#1]"))        return gaff::n;
        if (follows_rule(av, "[#7r6]([#6X3]=[#8])([#6X3]=[#8])[#6X4]")) return gaff::n;
        if (follows_rule(av, "[#7r6]([#6X3]=[#8])([#6X3]=[#8])[#1]"))   return gaff::n;
        if (follows_rule(av, "[#7r5]([#6X3]=[#8])([#6X3])[#6X4]"))      return gaff::n;
        if (follows_rule(av, "[#7r5]([#6X3]=[#8])([#6X3])[#1]"))        return gaff::n;
        if (follows_rule(av, "[#7]-[CX3]=[OA]"))                        return gaff::n;
        if (follows_rule(av, "[#7]-[CX3]=[SA]"))                        return gaff::n;

        // Amines
        if (follows_rule(av, "[#7H2]-[*r6;a]"))    return gaff::nh;
        if (follows_rule(av, "[NA]-[*R;a]"))       return gaff::nh;
        if (follows_rule(av, "[#7H2]-[*R;a]"))     return gaff::nh;

        // SP3 Nitrogen, strong evidence
        if (follows_rule(av, "[#7r6]([#6])([#6])[#1]"))   return gaff::n3;

        if (av.is_aromatic()) return gaff::na;

        if (follows_rule(av, "[#7r5]([#6])([#6])[#1]"))   return gaff::n3;
        if (follows_rule(av, "[#7r5]([#6])([#6])[#6]"))   return gaff::n3;
        if (follows_rule(av, "[#7r6]([#6])([#6])[#1]"))   return gaff::n3;
        if (follows_rule(av, "[#7r6]([#6])([#6])[#6]"))   return gaff::n3;

        if (follows_rule(av, "[#7;R](-*)=*-*"))   return gaff::na;

        return gaff::n3;
    case 4:
        return gaff::n4;
    default: // no idea
        return gaff::n4;
    }
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
        if (sssrs_count != 0 && av.is_aromatic()) return gaff::pb;
        if (sssrs_count != 0) {
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
        if (av.total_hydrogens() != 0)     return gaff::sh;
        if (follows_rule(av, "[SX2]=*"))   return gaff::s2;
        if (follows_rule(av, "[SX2]#*"))   return gaff::s2;
        return gaff::ss;
    case 3:
        if (follows_rule(av, "[SX3](=*)-*=*"))   return gaff::sx;
        if (follows_rule(av, "[SX3](=*)-*#*"))   return gaff::sx;
        return gaff::s4;
    case 4: // check for conjugation
        if (follows_rule(av, "[SX4](=*)-*=*"))   return gaff::sy;
        if (follows_rule(av, "[SX4](=*)-*#*"))   return gaff::sy;
        // fall through
    case 5:
    case 6:
    default: // no idea
        return gaff::s6;
    }
}

static size_t type_atom(const AtomVertex& av) {

    auto res = av.br().topology().residue_for_atom(av);
    if (res) {
        auto std_res = standard_residues.find(res->name());
        if (std_res != standard_residues.end()) {
            auto std_res_atom = std_res->second.find(av.name());
            if (std_res_atom != std_res->second.end()) {
                return atomtype_id_for_name<GAFF>(std_res_atom->second);
            }
        }
    }

    switch(av.atomic_number()) {
    case Element::H:
        return type_hydrogen(av);
        break;
    case Element::C:
        return type_carbon(av);
        break;
    case Element::N:
        return type_nitrogen(av);
        break;
    case Element::O:
        return type_oxygen(av);
        break;
    case Element::P:
        return type_phosphorus(av);
        break;
    case Element::S:
        return type_sulfur(av);
        break;
    case Element::F:
        return gaff::f;
        break;
    case Element::Cl:
        return gaff::cl;
        break;
    case Element::Br:
        return gaff::br;
        break;
    case Element::I:
        return gaff::i;
        break;
    default:
        return gaff::X;
        break;
    }
}

size_t GAFF::add_atom(size_t new_idx) {
    return add_atom_(new_idx);
}

size_t GAFF::add_atom_(size_t new_idx) {
    auto av = mol_[new_idx];

    if (new_idx >= this->size()) {
        this->resize(new_idx + 1);
    }

    set(new_idx, type_atom(av));

    return get(new_idx);
}

void GAFF::add_bond(size_t idx1, size_t idx2, Bond::Order /*bo*/) {
    set(idx1, type_atom(mol_[idx1]));
    set(idx2, type_atom(mol_[idx2]));
}

template<> std::string Spear::atomtype_name_for_id<GAFF>(size_t id) {
    return gaff_unmask[id];
}

template<> size_t Spear::atomtype_id_for_name<GAFF>(const std::string& name) {
    // I am lazy
    for (size_t i = 0; i < atomtype_id_count<GAFF>(); ++i) {
        if (atomtype_name_for_id<GAFF>(i) == name) {
            return i;
        }
    }
	return 0;
}

template<> size_t Spear::atomtype_id_count<GAFF>() {
	return gaff::Lr + 1;
}

template<> double Spear::van_der_waals<GAFF>(size_t /*id*/) {
	return 0.0;
}
