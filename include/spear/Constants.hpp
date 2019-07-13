#ifndef SPEAR_CONSTANTS_HPP
#define SPEAR_CONSTANTS_HPP

#include <vector>
#include <map>
#include <string>
#include <cstdint>

namespace Spear {
namespace Element {
enum Symbol : uint64_t {
    LP = 0,
    H = 1, D = 1, T = 1,                                                He,
    Li, Be,                                           B,  C,  N, O, F,  Ne,
    Na, Mg,                                           Al, Si, P, S, Cl, Ar,
    K,  Ca, Sc, Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr, Y,  Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
    Cs, Ba, La,
                Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu,
                Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac,
                Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
                Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Ml, Fc, Lv, Ts, Og,
};

const std::vector<std::string> Name {
    "LP",
    "H",                                                                                                  "He",
    "Li", "Be",                                                             "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg",                                                             "Al", "Si", "P",  "S",  "Cl", "Ar",
    "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs", "Ba", "La",
                      "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                      "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr",  "Ra", "Ac", 
                      "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                      "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Ml", "Fc", "Lv", "Ts", "Og",
};

// Taken from https://en.wikipedia.org/wiki/Covalent_radius
const std::vector<double> CovalentRadius {
    0.00,
    0.31,                                                                                                 0.28,
    1.28, 0.96,                                                             0.84, 0.68, 0.71, 0.66, 0.57, 0.58,
    1.66, 1.41,                                                             1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
    // The rest (below) need to be updated
    1.33, 0.99, 1.44, 1.47, 1.33, 1.35, 1.35, 1.34, 1.33, 1.50, 1.52, 1.45, 1.22, 1.17, 1.21, 1.22, 1.21, 1.16,
    1.47, 1.12, 1.78, 1.56, 1.48, 1.47, 1.35, 1.40, 1.45, 1.50, 1.59, 1.69, 1.63, 1.46, 1.46, 1.47, 1.40, 1.40,
    1.67, 1.34, 1.87,
                      1.83, 1.82, 1.81, 1.80, 1.80, 1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72,
                      1.57, 1.43, 1.37, 1.35, 1.37, 1.32, 1.50, 1.50, 1.70, 1.55, 1.54, 1.54, 1.68, 0.00, 1.50,
    0.00, 1.90, 1.88,
                      1.79, 1.61, 1.58, 1.55, 1.53, 1.51, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
                      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
};

namespace {
    std::map<std::string, Symbol> make_symbol_map_() {
        std::map<std::string, Symbol> map;
        for (size_t i = LP; i <= Og; ++i) {
            map.insert({Name[i], static_cast<Symbol>(i)});
        }
        return map;
    }
}

const auto SymbolForName = make_symbol_map_();

}

namespace Bond {
/// Keep in sync with Chemfiles::Connectivity::BondOrder
enum Order : uint64_t {
    UNKNOWN = 0,
    SINGLE = 1,
    DOUBLE = 2,
    TRIPLE = 3,
    QUADRUPLE = 4,
    QINTUPLET = 5,

    // space for more bond types if needed
    AMIDE = 254,
    AROMATIC = 255,
};
}

namespace Atom {
enum Chirality : uint64_t {
    UNSPECIFIED = 0,
    NONE,
    CIP_R,
    CIP_S,
    OTHER,
};
}

}

#endif
