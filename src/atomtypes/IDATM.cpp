// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include "spear/atomtypes/IDATM.hpp"
#include <map>
#include <string>
#include <locale>
#include <unordered_set>

//FIXME
#include <iostream>

#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

#include "spear/Geometry.hpp"

using namespace Spear;

#define RNA_CORE {"P","Pac"},{"OP1","O3-"},{"OP2","O3-"},\
                 {"O5'","O3"},{"C5'","C3"},{"C4'","C3"},{"O4'","O3"},{"C3'","C3"},{"O3'","O3"},{"C2'","C3"},{"O2'","O3"},{"C1'","C3"},{"N1","Npl"}
#define AMINO_ACID_CORE {"OXT","O2-"},{"CA","C3"},{"C","C2"},{"O","O2"},{"N","Npl"}

const std::map<std::string, std::map<std::string,std::string>> standard_residues {
    {"A",{RNA_CORE,{"N3","Npl"},{"N6","Npl"},{"N7","Npl"},{"N9","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"C8","Car"}}},
    {"C",{RNA_CORE,{"N3","Npl"},{"N4","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"O2","O2"}}},
    {"G",{RNA_CORE,{"N2","Npl"},{"N3","Npl"},{"N7","Npl"},{"N9","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"C8","Car"},{"O6","O2"}}},
    {"U",{RNA_CORE,{"N3","Npl"},{"C2","Car"},{"C4","Car"},{"C5","Car"},{"C6","Car"},{"O2","O2"},{"O4","O2"}}},
    {"ALA",{AMINO_ACID_CORE,{"CB","C3"}}},
    {"ARG",{AMINO_ACID_CORE,{"CB","C3"},{"CD","C3"},{"CG","C3"},{"CZ","C2"},{"NE","Ng+"},{"NH1","Ng+"},{"NH2","Ng+"}}},
    {"ASN",{AMINO_ACID_CORE,{"CB","C3"},{"CG","C2"},{"ND2","Npl"},{"OD1","O2"}}},
    {"ASP",{AMINO_ACID_CORE,{"CB","C3"},{"CG","Cac"},{"OD1","O2-"},{"OD2","O2-"}}},
    {"CYS",{AMINO_ACID_CORE,{"CB","C3"},{"SG","S3"}}},
    {"GLN",{AMINO_ACID_CORE,{"CB","C3"},{"CD","C2"},{"CG","C3"},{"NE2","Npl"},{"OE1","O2"}}},
    {"GLU",{AMINO_ACID_CORE,{"CB","C3"},{"CD","Cac"},{"CG","C3"},{"OE1","O2-"},{"OE2","O2-"}}},
    {"GLY",{AMINO_ACID_CORE}},
    {"HIS",{AMINO_ACID_CORE,{"CB","C3"},{"CD2","Car"},{"CE1","Car"},{"CG","Car"},{"ND1","Npl"},{"NE2","N2"}}},
    {"ILE",{AMINO_ACID_CORE,{"CB","C3"},{"CD1","C3"},{"CG1","C3"},{"CG2","C3"},}},
    {"LEU",{AMINO_ACID_CORE,{"CB","C3"},{"CD1","C3"},{"CD2","C3"},{"CG","C3"},}},
    {"LYS",{AMINO_ACID_CORE,{"CB","C3"},{"CD","C3"},{"CE","C3"},{"CG","C3"},{"NZ","N3+"},}},
    {"MET",{AMINO_ACID_CORE,{"CB","C3"},{"CE","C3"},{"CG","C3"},{"SD","S3"}}},
    {"MSE",{AMINO_ACID_CORE,{"CB","C3"},{"CE","C3"},{"CG","C3"},{"SE","Se"}}},
    {"PHE",{AMINO_ACID_CORE,{"CB","C3"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"}}},
    {"PRO",{AMINO_ACID_CORE,{"CB","C3"},{"CD","C3"},{"CG","C3"}}},
    {"SER",{AMINO_ACID_CORE,{"CB","C3"},{"OG","O3"}}},
    {"THR",{AMINO_ACID_CORE,{"CB","C3"},{"CG2","C3"},{"OG1","O3"}}},
    {"TRP",{AMINO_ACID_CORE,{"CB","C3"},{"CD1","Car"},{"CD2","Car"},{"CE2","Car"},{"CE3","Car"},{"CG","Car"},{"CH2","Car"},{"CZ2","Car"},{"CZ3","Car"},{"NE1","Npl"}}},
    {"TYR",{AMINO_ACID_CORE,{"CB","C3"},{"CD1","Car"},{"CD2","Car"},{"CE1","Car"},{"CE2","Car"},{"CG","Car"},{"CZ","Car"},{"OH","O3"}}},
    {"VAL",{AMINO_ACID_CORE,{"CB","C3"},{"CG1","C3"},{"CG2","C3"}}},
    {"HOH",{{"O","O3"}, {"OW", "O3"}, {"WO", "O3"}}},
};

// angle values used to discriminate between hybridization states
constexpr double angle23val1 = 115.0;
constexpr double angle23val2 = 122.0;
constexpr double angle12val = 160.0;

// bond length cutoffs from hybridization discrimination
// p3... = pass 3 cutoffs; p4... = pass 4 cutoffs
constexpr double p3c1c1 = 1.22;
constexpr double p3c2c = 1.41;
constexpr double p3c2n = 1.37;
constexpr double p3n1c1 = 1.20;
constexpr double p3n3c = 1.38;
constexpr double p3n3n3 = 1.43;
constexpr double p3n3n2 = 1.41;
constexpr double p3o2c2 = 1.30;
constexpr double p3o2as = 1.685;
constexpr double p3s2c2 = 1.76;
constexpr double p3s2as = 2.11;
constexpr double p4c3c = 1.53;
constexpr double p4c3n = 1.46;
constexpr double p4c3o = 1.44;
constexpr double p4n2c = 1.38;
constexpr double p4n2n = 1.32;
constexpr double p4c2c = 1.42;
constexpr double p4c2n = 1.41;
constexpr double p4ccnd = 1.45;
constexpr double p8cc2n2h = 1.367;
constexpr double p8nn2n2h = 1.326;
constexpr double p8cn2n2h = 1.367;

enum idatm : size_t {
    Ac, Ag, Al, Am, Ar, As, At, Au, B, Ba, Be, Bh, Bi, Bk, Br, C, C1, C1m, C2,
    C3, Ca, Cac, Car, Cd, Ce, Cf, Cl, Cm, Co, Cr, Cs, Cu, D, Db, DC, Ds, Dy, Er,
    Es, Eu, F, Fe, Fm, Fr, Ga, Gd, Ge, H, HC, He, Hf, Hg, Ho, Hs, I, In, Ir, K,
    Kr, La, Li, Lr, Lu, Lw, Md, Mg, Mn, Mo, Mt, N, N1, N1p, N2, N2p, N3, N3p,
    Na, Nb, Nd, Ne, Ngp, Ni, No, Nox, Np, Npl, Ntr, O, O1, O1p, O2, O2m, O3,
    O3m, Oar, Oarp, Os, P, P3p, Pa, Pac, Pb, Pd, Pm, Po, Pox, Pr, Pt, Pu, Ra,
    Rb, Re, Rf, Rh, Rn, Ru, S,S2, S3, S3m, S3p, Sac, Sar, Sb, Sc, Se, Sg, Si,
    Sm, Sn, Son, Sr,Sxd, Ta, Tb, Tc, Te,Th, Ti, Tl, Tm, U,V, W, Xe, Y, Yb, Zn,
    Zr, unk
};

const std::unordered_map<std::string, size_t> idatm_mask{
    {"Ac", Ac},   {"Ag", Ag},   {"Al", Al},   {"Am", Am},  {"Ar", Ar},
    {"As", As},   {"At", At},   {"Au", Au},   {"B", B},    {"Ba", Ba},
    {"Be", Be},   {"Bh", Bh},   {"Bi", Bi},   {"Bk", Bk},  {"Br", Br},
    {"C", C},     {"C1", C1},   {"C1-", C1m}, {"C2", C2},  {"C3", C3},
    {"Ca", Ca},   {"Cac", Cac}, {"Car", Car}, {"Cd", Cd},  {"Ce", Ce},
    {"Cf", Cf},   {"Cl", C1},   {"Cm", Cm},   {"Co", Co},  {"Cr", Cr},
    {"Cs", Cs},   {"Cu", Cu},   {"D", D},     {"Db", Db},  {"DC", DC},
    {"Ds", Ds},   {"Dy", Dy},   {"Er", Er},   {"Es", Es},  {"Eu", Eu},
    {"F", F},     {"Fe", Fe},   {"Fm", Fm},   {"Fr", Fr},  {"Ga", Ga},
    {"Gd", Gd},   {"Ge", Ge},   {"H", H},     {"HC", HC},  {"He", He},
    {"Hf", Hf},   {"Hg", Hg},   {"Ho", Ho},   {"Hs", Hs},  {"I", I},
    {"In", In},   {"Ir", Ir},   {"K", K},     {"Kr", Kr},  {"La", La},
    {"Li", Li},   {"Lr", Lr},   {"Lu", Lu},   {"Lw", Lw},  {"Md", Md},
    {"Mg", Mg},   {"Mn", Mn},   {"Mo", Mo},   {"Mt", Mt},  {"N", N},
    {"N1", N1},   {"N1+", N1p}, {"N2", N2},   {"N2+", N2p},{"N3", N3},
    {"N3+", N3p}, {"Na", Na},   {"Nb", Nb},   {"Nd", Nd},  {"Ne", Ne},
    {"Ng+", Ngp}, {"Ni", Ni},   {"No", No},   {"Nox", Nox},{"Np", Np},
    {"Npl", Npl}, {"Ntr", Ntr}, {"O", O},     {"O1", O1},  {"O1+", O1p},
    {"O2", O2},   {"O2-", O2m}, {"O3", O3},   {"O3-", O3m},{"Oar", Oar},
    {"Oar+",Oarp},{"Os", Os},   {"P", P},     {"P3+", P3p},{"Pa", Pa},
    {"Pac", Pac}, {"Pb", Pb},   {"Pd", Pd},   {"Pm", Pm},  {"Po", Po},
    {"Pox", Pox}, {"Pr", Pr},   {"Pt", Pt},   {"Pu", Pu},  {"Ra", Ra},
    {"Rb", Rb},   {"Re", Re},   {"Rf", Rf},   {"Rh", Rh},  {"Rn", Rn},
    {"Ru", Ru},   {"S", S},     {"S2", S2},   {"S3", S3},  {"S3-", S3m},
    {"S3+", S3p}, {"Sac", Sac}, {"Sar", Sar}, {"Sb", Sb},  {"Sc", Sc},
    {"Se", Se},   {"Sg", Sg},   {"Si", Si},   {"Sm", Sm},  {"Sn", Sn},
    {"Son", Son}, {"Sr", Sr},   {"Sxd", Sxd}, {"Ta", Ta},  {"Tb", Tb},
    {"Tc", Tc},   {"Te", Te},   {"Th", Th},   {"Ti", Ti},  {"Tl", Tl},
    {"Tm", Tm},   {"U", U},     {"V", V},     {"W", W},    {"Xe", Xe},
    {"Y", Y},     {"Yb", Yb},   {"Zn", Zn},   {"Zr", Zr},  {"unk", unk},
};

const char* const idatm_unmask[]{
    "Ac",  "Ag",  "Al", "Am",  "Ar", "As",  "At",  "Au",   "B",   "Ba",  "Be",
    "Bh",  "Bi",  "Bk", "Br",  "C",  "C1",  "C1-", "C2",   "C3",  "Ca",  "Cac",
    "Car", "Cd",  "Ce", "Cf",  "Cl", "Cm",  "Co",  "Cr",   "Cs",  "Cu",  "D",
    "Db",  "DC",  "Ds", "Dy",  "Er", "Es",  "Eu",  "F",    "Fe",  "Fm",  "Fr",
    "Ga",  "Gd",  "Ge", "H",   "HC", "He",  "Hf",  "Hg",   "Ho",  "Hs",  "I",
    "In",  "Ir",  "K",  "Kr",  "La", "Li",  "Lr",  "Lu",   "Lw",  "Md",  "Mg",
    "Mn",  "Mo",  "Mt", "N",   "N1", "N1+", "N2",  "N2+",  "N3",  "N3+", "Na",
    "Nb",  "Nd",  "Ne", "Ng+", "Ni", "No",  "Nox", "Np",   "Npl", "Ntr", "O",
    "O1",  "O1+", "O2", "O2-", "O3", "O3-", "Oar", "Oar+", "Os",  "P",   "P3+",
    "Pa",  "Pac", "Pb", "Pd",  "Pm", "Po",  "Pox", "Pr",   "Pt",  "Pu",  "Ra",
    "Rb",  "Re",  "Rf", "Rh",  "Rn", "Ru",  "S",   "S2",   "S3",  "S3-", "S3+",
    "Sac", "Sar", "Sb", "Sc",  "Se", "Sg",  "Si",  "Sm",   "Sn",  "Son", "Sr",
    "Sxd", "Ta",  "Tb", "Tc",  "Te", "Th",  "Ti",  "Tl",   "Tm",  "U",   "V",
    "W",   "Xe",  "Y",  "Yb",  "Zn", "Zr",  "unk"
};

const std::unordered_set<size_t> c2impossibleTypes = {
    idatm::C3,  idatm::DC,  idatm::HC, idatm::N3,  idatm::N3p, idatm::O3,
    idatm::Cac, idatm::Pac, idatm::Sac, idatm::Son, idatm::C1,  idatm::S3,
    idatm::Nox, idatm::Sxd,
};

const std::unordered_set<size_t> possible_aromatic = {
    idatm::C2, idatm::Npl, idatm::S2,  idatm::O3, idatm::S3, idatm::N3,
    idatm::Oar, idatm::P, idatm::Car, idatm::N2, idatm::O2,
};

const std::unordered_set<size_t> sp3hybridized = {
    idatm::C3, idatm::N3p, idatm::N3, idatm::O3, idatm::S3p, idatm::S3,
    idatm::P3p
};

const std::unordered_set<size_t> oxygenConjugatedTetra = {
    idatm::N3p, idatm::N2p, idatm::Pac, idatm::Pox, idatm::Sac,
    idatm::Son, idatm::Sxd, idatm::Nox,
};

const std::map<Element::Symbol, std::map<Bond::Order, size_t>> bond_maps = {
    {Element::C, {{Bond::DOUBLE, idatm::C2}, {Bond::TRIPLE, idatm::C1},
                  {Bond::AROMATIC, idatm::Car}}},
    {Element::N, {{Bond::DOUBLE, idatm::N2}, {Bond::TRIPLE, idatm::N1},
                  {Bond::AROMATIC, idatm::Npl}, {Bond::AMIDE, idatm::Npl}}},
    {Element::O, {{Bond::DOUBLE, idatm::O2}, {Bond::TRIPLE, idatm::O1},
                  {Bond::AROMATIC, idatm::Oar}}},
    {Element::S, {{Bond::DOUBLE, idatm::S2}, {Bond::AROMATIC, idatm::Sar}}},
};

static bool check_dihedrals_for_planarity(const Molecule& mol, const std::set<size_t>& ring) {
    auto& dihedrals = mol.topology().dihedrals();
    // test dihedral angles
    for (const auto& d : dihedrals) {

        if (ring.count(d[0]) == 0 || ring.count(d[1]) == 0 ||
            ring.count(d[2]) == 0 || ring.count(d[3]) == 0) {
            continue;
        }

        auto dangle = dihedral(mol[d[0]].position(),
                               mol[d[1]].position(),
                               mol[d[2]].position(),
                               mol[d[3]].position());

        dangle *= 180 / 3.14159;

        // larger than 10 degrees, means NOT aromatic
        if (std::fabs(dangle) > 10 &&
            std::fabs(std::fabs(dangle) - 180) > 10) {
            return false;
        }
    }
    // aromatic
    return true;
}

static bool aromatic(const Molecule& mol, const std::set<size_t>& ring) {
    auto sum = 0.0;
    size_t bonds = 0;
    const std::vector<BondEdge> ring_bonds = mol.get_bonds_in(ring);
    for (const auto& bond : ring_bonds) {
        const auto atom1 = bond.source();
        const auto atom2 = bond.target();
        const auto e1 = atom1.atomic_number();
        const auto e2 = atom2.atomic_number();
        double d = distance(atom1.position(), atom2.position());
        if (e1 == 6 && e2 == 6) {
            bonds++;
            sum += (d - 1.397) * (d - 1.397);
        } else if ((e1 == 6 || e2 == 6) &&
                   (e1 == 7 || e2 == 7)) {
            bonds++;
            sum += (d - 1.338) * (d - 1.338);
        } else if (e1 == 7 && e2 == 7) {
            bonds++;
            sum += (d - 1.308) * (d - 1.308);
        } else if ((e1 == 6 || e2 == 6) &&
                   (e1 == 8 || e2 == 8)) {
            bonds++;
            sum += (d - 1.300) * (d - 1.300);
        }
    }
    if (bonds == 0) return false;
    double homas = 1.0 - (98.89 / static_cast<double>(bonds)) * sum;
    bool is_aromatic = true;
    if (homas < 0.271) {
        is_aromatic = check_dihedrals_for_planarity(mol, ring);
    }
    return is_aromatic;
}

static size_t freeOxygens(const Spear::AtomVertex& atom,
                          const std::vector<size_t>& heavys) {
    size_t freeOxygens = 0;
    for (auto neighbor : atom.neighbors()) {
        if (neighbor.atomic_number() == Element::O && heavys[neighbor] == 1) {
            ++freeOxygens;
        }
    }

    return freeOxygens;
}

static size_t assignBondOrderType(Element::Symbol atomic_number,
                                  Bond::Order bo,
                                  size_t current_type) {
    if (current_type == idatm::N2 && bo == Bond::DOUBLE) {
        return idatm::N1;
    }

    // Allenes 
    if (current_type == idatm::C2 && bo == Bond::DOUBLE) {
        return idatm::C1;
    }

    if (current_type != idatm::C3 && current_type != idatm::N3 &&
        current_type != idatm::O3) {
        return current_type;
    }

    auto number_map = bond_maps.find(atomic_number);
    if (number_map != bond_maps.end()) {
        auto order_map = number_map->second.find(bo);
        if (order_map != number_map->second.end()) {
            return order_map->second;
        }
    }

    return current_type;
}

IDATM::IDATM(const Molecule& mol, TypingMode mode) :
    mol_(mol),
    heavys_(mol.size(), 0),
    mapped_(mol.size(), false),
    redo(mol.size())
{
    resize(mol.size(), idatm::unk);
    switch (mode) {
    case GEOMETRY:
        name_ = "IDATM_geometry";
        type_atoms_3d_();
        break;
    case TOPOLOGY:
        name_ = "IDATM_topology";
        type_atoms_topo_();
        break;
    default:
        throw std::invalid_argument("Unknown mode");
        break;
    }
}

template<typename func>
void IDATM::apply_function_to_all_atoms(func&& f) {
    for(auto atom : mol_) {
        if (mapped_[atom]) {
            continue;
        }

         set(atom, f(this, atom));
    }
}

void IDATM::type_atoms_3d_() {
    infallible_();
    apply_function_to_all_atoms(std::mem_fn(&IDATM::valence_));
    apply_function_to_all_atoms(std::mem_fn(&IDATM::terminal_));
    apply_function_to_all_atoms(std::mem_fn(&IDATM::redo_));
    fix_C2_();
    charges_();
    aromatic_();
    fix_N2_();
    heteroaromatics_();
    apply_function_to_all_atoms(std::mem_fn(&IDATM::fix_tetrahedrals_));
    apply_function_to_all_atoms(std::mem_fn(&IDATM::fix_special_));
}

void IDATM::type_atoms_topo_() {

    infallible_();
    apply_function_to_all_atoms(std::mem_fn(&IDATM::valence_topo_));

    size_t aromatic_count = 0;
    for (auto bond : mol_.bonds()) {
        if (bond.order() == Bond::AROMATIC) {
            ++aromatic_count;

            // Spoof the aromatic ring list
            aromatic_ring_sizes_.insert({bond.source(), 0});
            aromatic_ring_sizes_.insert({bond.target(), 0});
        }

        set(bond.source(),
            assignBondOrderType(bond.source().atomic_number(),
                                bond.order(),
                                get(bond.source())
            )
        );

        set(bond.target(),
            assignBondOrderType(bond.target().atomic_number(),
                                bond.order(),
                                get(bond.target())
            )
        );
    }

    fix_C2_();
    charges_();

    if (aromatic_count < 5 && dimensionality(mol_.positions()) == 3) {
        aromatic_ring_sizes_.clear();
        aromatic_();
    }

    for (auto bond : mol_.bonds()) {
        if (get(bond.source()) == idatm::N3 && get(bond.target()) == idatm::Car) {
            set(bond.source(), idatm::Npl);
        }

        if (get(bond.target()) == idatm::N3 && get(bond.source()) == idatm::Car){
            set(bond.target(), idatm::Npl);
        }
    }

    fix_N2_();
    heteroaromatics_();
    apply_function_to_all_atoms(std::mem_fn(&IDATM::fix_tetrahedrals_));
    apply_function_to_all_atoms(std::mem_fn(&IDATM::fix_special_));
}

void IDATM::infallible_() {
    for (const auto atom : mol_) {
        if (atom.atomic_number() == 1) { // Hydrogen and deuterium
            bool bondedToCarbon = false;
            for (auto neighbor : atom.neighbors()) {
                if (neighbor.atomic_number() == 6) {
                    bondedToCarbon = true;
                    break;
                }
            }

            // TODO: Update once Chemfiles is fixed
            bool isHyd = mol_.topology()[atom].mass() != 2.0;

            set(atom, bondedToCarbon ? (isHyd ? HC: DC)
                                     : (isHyd ? H : D)
            );

            // Perfect evidence
            mapped_[atom] = true;

            continue;
        }

        auto neighs = atom.neighbors();

        auto heavys = std::count_if(neighs.begin(), neighs.end(),
            [](AtomVertex n) {return n.atomic_number() > 1;}
        );

        heavys_[atom] = static_cast<std::size_t>(heavys);
    }

    // Use templates for "infallible" typing of standard residues
    for (auto& residue : mol_.topology().residues()) {
        auto it = standard_residues.find(residue.name());
        if (it != standard_residues.end()) {

            for (auto i : residue) {
                // Hydrogens are mapped
                if (mapped_[i]) {
                    continue;
                }

                mapped_[i] = true;
                auto a = mol_[i];

                // is it the N-terminal residue ?
                if (a.name() == "N" && heavys_[i] == 1) {
                    set(i, idatm::N3p);
                    continue;
                }

                // is it C-terminal ?
                if (a.name() == "C" && freeOxygens(a, heavys_) == 2) {
                    set(i, idatm::Cac);
                    continue;
                }

                auto it2 = it->second.find(a.name());
                if (it2 == it->second.end()) {
                    /// FIXME Using logging instead
                    std::cerr << "Cannot find atom name " + a.name()
                              << " of template residue " + residue.name();
                    mapped_[i] = false;
                    continue;
                }
                set(i, idatm_mask.at(it2->second));
            }
        }
    }

    for (auto atom : mol_) {
        // Hydrogens, standard residues
        if (mapped_[atom]) {
            continue;
        }

        // undifferentiated types
        auto element = atom.atomic_number();
        if ((element >= 2 && element <= 4) ||
            (element >= 10 && element <= 14) ||
            element >= 17) {
            set(atom, idatm_mask.at(Element::Name[atom.atomic_number()]));
            mapped_[atom] = true; // infallible type
            continue;
        }
    }
}

size_t IDATM::valence_(const AtomVertex& atom) {

    auto element = atom.atomic_number();
    auto valence = atom.degree();
    auto freeOs = freeOxygens(atom, heavys_);

    if (valence == 4) { // assume tetrahedral
        if (element == Element::C) {
            return idatm::C3; // must be sp3 carbon
        } else if (element == Element::N) {
            return freeOs >= 1 ? idatm::Nox
                               : idatm::N3p;
        } else if (element == Element::P) {
            if (freeOs >= 2) { // phostphate
                return idatm::Pac;
            } else if (freeOs == 1) { // P-oxide
                return idatm::Pox;
            } else { // Formally positive SP3 phospohrus
                return idatm::P3p;
            }
        } else if (element == Element::S) {
            if (freeOs >= 3) { // Sulfate
                return idatm::Sac;
            } else if (freeOs >= 1) { // Sulfone
                return idatm::Son;
            } else { // We don't know! Other!
                return idatm::S;
            }
        }
    } else if (valence == 3) {
        auto avgAngle = 0.0;
        for (size_t n1 = 0; n1 < 3; ++n1) {
            for (size_t n2 = n1 + 1; n2 < 3; ++n2) {
                avgAngle += angle(atom[n1].position(),
                                  atom.position(),
                                  atom[n2].position());
            }
        }
        avgAngle /= 3.0;
        avgAngle *= 180 / 3.14149;

        if (element == Element::C) {
            if (avgAngle < angle23val1) { // angle significantly < 120?
                return idatm::C3; // Then tetrahedral
            } else { // Most likely trigonal planar (some expceptions)
                return freeOs >= 2 ? idatm::Cac
                                   : idatm::C2;
            }
        } else if (element == Element::N) {
            if (avgAngle < angle23val1) {
                return idatm::N3; // likely tetrahedral
            } else { 
                return freeOs >= 2 ? idatm::Ntr
                                                : idatm::Npl;
            }
        } else if (element == Element::S) { // sulfoxide or formally positive sp3 S
            return freeOs > 0 ? idatm::Sxd
                              : idatm::S3p;
        }
    } else if (valence == 2) {
        double ang = angle(atom[0].position(),
                           atom.position(),
                           atom[1].position()) * 180 / 3.14159;

        if (element == Element::C) {
            if (ang < angle23val1) { // could be tetralhedral, let's redo
                redo[atom] = 1;
                return idatm::C3;
            } else if (ang < angle12val) { // Signficantly less than 180o
                if (ang < angle23val2) { // Angle is too small for comfort
                    redo[atom] = 3;
                }
                return idatm::C2;
            } else { // It's near 180o, so sp carbon
                return idatm::C1;
            }
        } else if (element == Element::N) {
            if (ang < angle23val1) { // Check if its tetralhedral, redo
                return idatm::N3;
                redo[atom] = 2;
            } else { // Is it near 180o?
                return ang < angle12val ? idatm::Npl
                                        : idatm::N1;
            }
        } else if (element == Element::O) { // Valence of two = tetrahedral
            return idatm::O3;
        } else if (element == Element::S) {
            return idatm::S3;
        }
    }

    // For valence 1, ensure that a default type gets assigned.
    return idatm_mask.at(Element::Name[atom.atomic_number()]);
}

size_t IDATM::valence_topo_(const AtomVertex& atom) {
    auto freeOs = freeOxygens(atom, heavys_);
    auto valence = atom.degree();

    switch (atom.atomic_number()) { 
    case Element::C:
        return idatm::C3;
    case Element::N:
        if (valence == 4) {
            return freeOs >= 1 ? idatm::Nox
                               : idatm::N3p;
        }
        if (valence == 3) {
            return freeOs >= 2 ? idatm::Ntr
                               : idatm::N3;
        }
        return idatm::N3;
    case Element::O:
        return idatm::O3;
    case Element::P:
        if (valence == 4) {
            if (freeOs >= 2) { // phosphate
                return idatm::Pac;
            }
            if (freeOs == 1) { // P-oxide
                return idatm::Pox;
            }
            // Formally positive SP3 phosphorus
            return idatm::P3p;
        }
        return idatm::P;
    case Element::S:
        if (valence == 4) {
            if (freeOs >= 3) { // Sulfate
                return idatm::Sac;
            }
    
            if (freeOs >= 1) { // Sulfone
                return idatm::Son;
            }

            // We don't know! Other!
            return idatm::S;
        }
        if (valence == 3) {
            return freeOs > 0 ? idatm::Sxd
                              : idatm::S3p;
        }
        return idatm::S3;
    default:
        break;
    }

    return get(atom);
}

size_t IDATM::terminal_(const AtomVertex& atom) {
    if (atom.degree() != 1) {
        return get(atom);
    }

    auto bondee = atom[0];
    auto len = distance(atom.position(), bondee.position());
    auto bondeeType = get(bondee);

    switch(get(atom)) {
    case idatm::C: // Default Carbon
        if (len <= p3c1c1 && bondeeType == idatm::C1) { // Check length
            return idatm::C1;
        }
        if (len <= p3c2c && bondee.atomic_number() == Element::C) {
            return idatm::C2;
        }
        if (len <= p3c2n && bondee.atomic_number() == Element::N) {
            return idatm::C2;
        }
        return idatm::C3;
    case idatm::N: // Default Nitrogen
        if (len <= p3n1c1 && bondeeType == idatm::C1) {
            return idatm::N1;
        }
        if (len > p3n3c &&
            (bondeeType == idatm::C2 ||
             bondeeType == idatm::C3)) {
            return idatm::N3;
        }
        if ((len > p3n3n3 && bondeeType == idatm::N3) ||
            (len > p3n3n2 && bondeeType == idatm::Npl)) {
            return idatm::N3;
        }
        return idatm::Npl;
    case idatm::O: // Default Oxygen
        if (bondeeType == idatm::Cac ||
            bondeeType == idatm::Pac ||
            bondeeType == idatm::Sac ||
            bondeeType == idatm::Ntr) {
            return idatm::O2m;
        }
        if (bondeeType == idatm::Nox ||
            bondeeType == idatm::Pox ||
            bondeeType == idatm::Son ||
            bondeeType == idatm::Sxd) {
            return idatm::O2;
        }
        if (len <= p3o2c2 && bondee.atomic_number() == Element::C) {
            if (!mapped_[bondee])
                set(bondee, idatm::C2);
            redo[bondee] = 0; // Don't redo bondee, we are sure.
            return idatm::O2;
        }
        if (len <= p3o2as && bondee.atomic_number() == Element::As) {
            return idatm::O2;
        }
        return idatm::O3;
    case idatm::S: // Default Sulfur
        if (bondee.atomic_number() == Element::P) {
            return idatm::S2;
        }
        if (len <= p3s2c2 && bondee.atomic_number() == Element::C) {
            if (!mapped_[bondee])
                set(bondee, idatm::C2);
            redo[bondee] = 0; // Don't redo bondee, we are sure.
            return idatm::S2;
        }
        if (len <= p3s2as && bondee.atomic_number() == Element::As) {
            return idatm::S2;
        }
        return idatm::S3;
    default:
        return get(atom);
    }
}

size_t IDATM::redo_(const AtomVertex& atom) {
    if (redo[atom] == 0) return get(atom);

    bool c3able = false;
    for (auto bondee : atom.neighbors()) {
        auto len = distance(atom.position(), bondee.position());
        auto bondeeElement = bondee.atomic_number();

        if (redo[atom] == 1) { // Tetrahedral or planar carbon?
            if ((len > p4c3c && bondeeElement == Element::C) ||
                (len > p4c3n && bondeeElement == Element::N) ||
                (len > p4c3o && bondeeElement == Element::O)) {
                return idatm::C3; // Done for this redo
            }
            if ((len <= p4c2c && bondeeElement == Element::C) ||
                (len <= p4c2n && bondeeElement == Element::N)) {
                return idatm::C2;
            }
        } else if (redo[atom] == 2) { // Tetrahedral or planar nitrogen?
            if ((len <= p4n2c && bondeeElement == Element::C) ||
                (len <= p4n2n && bondeeElement == Element::N)) {
                return idatm::Npl; // Done for this redo
            }
        } else {
            if ((len <= p4c2c && bondeeElement == Element::C) ||
                (len <= p4c2n && bondeeElement == Element::N)) {
                return idatm::C2;
            }
            if ((len > p4c3c && bondeeElement == Element::C) ||
                (len > p4c3n && bondeeElement == Element::N) ||
                (len > p4c3o && bondeeElement == Element::O)) {
                c3able = true;
            }

            if (len > p4ccnd && bondeeElement == Element::C) c3able = true;
        }
    }

    if (c3able) return idatm::C3;
    return get(atom);
}

void IDATM::fix_C2_() {
    for (auto atom : mol_) {
        if (get(atom) != idatm::C2) continue;

        if (freeOxygens(atom, heavys_) >= 2) {
            set(atom, idatm::Cac);
            continue;
        }

        if (mapped_[atom]) continue;

        bool c2possible = false;
        for (auto bondee : atom.neighbors()) {
            if (c2impossibleTypes.count(get(bondee)) == 0) {
                c2possible = true;
                break;
            }
        }
        if (!c2possible) set(atom, idatm::C3);
    }
}

void IDATM::charges_() {
    for (auto atom : mol_) {
        if (get(atom) == idatm::N3) {
            if (mapped_[atom]) continue;
            bool positive = true;
            auto c2_index = mol_.size();
            for (auto bondee : atom.neighbors()) {
                if (get(bondee) != idatm::C3 &&
                    get(bondee) != idatm::H &&
                    get(bondee) != idatm::D) {
                    positive = false;
                }
                if (get(bondee) == idatm::C2) {
                    c2_index = bondee;
                }
            }
            if (positive) return set(atom, idatm::N3p);

            if (c2_index != mol_.size() && freeOxygens(mol_[c2_index], heavys_) == 1) {
                return set(atom, idatm::Npl);
            }
        }

        if (get(atom) == idatm::C2) {
            int numNpls = 0;
            bool hasN2 = false;
            for (auto bondee : atom.neighbors()) {
                if ((get(bondee) == idatm::Npl && !mapped_[bondee]) ||
                    (get(bondee) == idatm::N3  && heavys_[bondee] == 1 && !mapped_[bondee]) ||
                    get(bondee) == idatm::Ngp) {
                    // Ng+ possible through template typing
                    numNpls++;
                }
                if (get(bondee) == idatm::N2 && heavys_[bondee] <= 2 && !mapped_[bondee]) {
                    hasN2 = true;
                }
            }

            bool noplus = false;
            if (numNpls == 3 || (numNpls == 2 && hasN2)) {
                for (auto bondee : atom.neighbors()) {
                    //if (atom_types_[bondee] != idatm::Npl) continue;
                    if (mapped_[bondee]) continue;

                    set(bondee, idatm::Ngp);
                    for (auto bondee2 : mol_[bondee].neighbors()) {
                        if ((get(bondee2) == idatm::C2 || get(bondee2) == idatm::Npl) &&
                            bondee2 != atom) {
                            set(bondee, idatm::Npl);
                            noplus = true;
                            break;
                        }
                    }
                }
            }
            // Reset!
            if (noplus) {
                for (auto bondee : atom.neighbors()) {
                    if (mapped_[bondee]) continue;
                    if (get(bondee) == idatm::Ngp) {
                        set(bondee, idatm::Npl);
                    }
                }
            }
        }

        if (get(atom) == idatm::Cac) {
            for (auto bondee : atom.neighbors()) {
                if (mapped_[bondee]) continue;
                if (bondee.atomic_number() == 8 && heavys_[bondee] == 1) {
                    set(bondee, idatm::O2m);
                }
            }
        }
    }
}

void IDATM::aromatic_() {
    auto rings = mol_.smallest_set_of_smallest_rings();
    for (const auto& ring : rings) {
        bool planarTypes = true;
        size_t c3_count = 0;
        for (auto atom : ring) {
            if (possible_aromatic.count(get(atom)) == 0) {
                if (get(atom) == idatm::C3 && mol_[atom].degree() == 2) {
                    ++c3_count;
                    continue;
                }
                // O3/S3/N3/C3 will be changed to Oar/S2/Npl/Car if
                // ring is found to be aromatic
                planarTypes = false;
                break;
            }

            if ((get(atom) == idatm::C2 ||
                 get(atom) == idatm::C3 ||
                 get(atom) == idatm::O3 ||
                 get(atom) == idatm::S3 ||
                 get(atom) == idatm::N3) &&
                mapped_[atom]) {
                planarTypes = false;
                break;
            }
        }
        if (!planarTypes) continue;

        // Special 0T8 fix
        if (ring.size() == 5 && c3_count == 3) {
            bool has_npl = false, has_c2 = false;
            for (auto atom : ring) {
                if (get(atom) == idatm::C2)  has_c2 = true;
                if (get(atom) == idatm::Npl) has_npl = true;
            }
            if (has_c2 && has_npl) continue;
        }

        if (!aromatic(mol_, ring)) {
            continue;
        }

        // Correct types to be aromatic
        for (auto atom : ring) {
            aromatic_ring_sizes_.insert({atom, ring.size()});
            switch(get(atom)) {
            case idatm::C2: set(atom, idatm::Car); break;
            case idatm::C3: set(atom, idatm::Car); break;
            case idatm::O3: set(atom, idatm::Oar); break;
            case idatm::O2: set(atom, idatm::Oarp); break;
            case idatm::S3: set(atom, idatm::S2); break;
            case idatm::S2: set(atom, idatm::Sar); break;
            case idatm::N3: set(atom, idatm::Npl); break;
            default: break;
            }
        }
    }
}

void IDATM::fix_N2_() {
    for (auto atom : mol_) {
        if (mapped_[atom]) continue;
        if (get(atom) != idatm::Npl) continue;
        if (heavys_[atom] != 2) continue;
        if (atom.degree() > 2) continue;

        // are both bonded heavy atoms sp3?  If so -> Npl
        bool bothSP3 = true;
        bool bothC = true;
        bool bothCar = true;
        bool bothN = true;
        bool other = false;
        for (auto bondee : atom.neighbors()) {
            if (sp3hybridized.count(get(bondee)) == 0) {
                bothSP3 = false;
                // don't break out, need to check both
                // bonded heavies' element types
            }

            if (bondee.atomic_number() == 6) {
                bothN = false;
                if (get(bondee) != idatm::Car) {
                    bothCar = false;
                }
            } else if (bondee.atomic_number() == 7) {
                bothC = false;
                bothCar = false;
            } else {
                other = true;
                bothN = false;
                bothC = false;
                bothCar = false;
            }
        }
        if (bothSP3) continue;

        // Bridging nitrogen between two aromatic rings!
        if (bothCar && aromatic_ring_sizes_.find(atom) == aromatic_ring_sizes_.end()) {
            continue;
        }

        if (other) {
            // Metal bond most likely. (C(=O)N:---M), keep Npl
            continue;
        }
        auto avgLen = 0.0;
        for (auto bondee : atom.neighbors()) {
            auto len = distance(atom.position(), bondee.position());
            avgLen += len;
        }
        avgLen /= 2.0;
        if (bothC) {
            if (avgLen <= p8cc2n2h) set(atom, idatm::N2);
        } else if (bothN) {
            if (avgLen <= p8nn2n2h) set(atom, idatm::N2);
        } else if (avgLen <= p8cn2n2h) {
            // one N, one C
            set(atom, idatm::N2);
        }
    }
}

void IDATM::heteroaromatics_() {
    // Do two passes
    // First, assign Oar+ and Sar
    for (auto pair : aromatic_ring_sizes_) {
        if (mapped_[pair.first]) continue;

        // aromatic sulfur, always assign??
        if (mol_[pair.first].atomic_number() == 16) {
            set(pair.first, idatm::Sar);
        }

        // Formally positive oxygen
        if (get(pair.first) == idatm::Oar && pair.second % 2 == 0) {
            set(pair.first, idatm::Oarp);
        }
    }

    // Now we need to retype only Npl
    for (auto pair : aromatic_ring_sizes_) {
        if (mapped_[pair.first]) continue;

        if (get(pair.first) != idatm::Npl &&
            get(pair.first) != idatm::N2) continue;

        // change 2-substituted Npl in 6-membered aromatic ring to N2
        // and a 3+-substituted Npl in 6-membered aromatic ring to N2+
        if (pair.second == 6) {
            set(pair.first, heavys_[pair.first] == 2 ? idatm::N2
                                                     : idatm::N2p
            );
            continue;
        }

        auto bound_to_car = false;
        auto bound_to_npl = false;
        auto bound_to_oar = false;
        auto bound_to_sar = false; // Note: changed from just S?

        for (auto bondee : mol_[pair.first].neighbors()) {
            switch (get(bondee)) {
            case idatm::Npl: bound_to_npl = true; break;
            case idatm::Car: bound_to_car = true; break;
            case idatm::Oar: bound_to_oar = true; break;
            case idatm::Oarp: bound_to_oar= true; break;
            case idatm::Sar: bound_to_sar = true; break;
            default: 
                break;
            }
        }

        // change aromatic Npl bound to Npl and Car to N2
        if (bound_to_npl && bound_to_car) {
            set(pair.first, idatm::N2);
        }

        // change 2/3-substituted Npl bound to S or O in 5-ring to N2
        if (pair.second == 5 && (bound_to_oar || bound_to_sar) && 
            (heavys_[pair.first] == 2 || heavys_[pair.first] == 3)) {
            set(pair.first, idatm::N2);
        }

        // This is odd, we change the first Npl in 5-ring
        // Npl(-H)-Car-Npl(-C,-C,-C) to N2
        // Similar to HIS, but this does not match intuitive atome types...
        if (pair.second == 5 && bound_to_car &&
            (mol_[pair.first].degree() > heavys_[pair.first])) {
            set(pair.first, idatm::N2);
        }
    }
}

size_t IDATM::fix_tetrahedrals_(const AtomVertex& atom) {
    if (heavys_[atom] != 1) return get(atom);

    auto bondee = atom[0];

    // resonance equivalent terminal oxygen on tetrahedral
    // center (phosphate, sulfate, sulfone...)
    if (atom.atomic_number() == Element::O &&
        oxygenConjugatedTetra.count(get(bondee)) != 0) {
        return idatm::O3m;
    }

    // change terminal C3 bound to C1, which in turn is
    // bound to any C but C1 (Ligand ID: 1DJ, Atom Name:
    // CAA) to C1
    if (get(atom) == idatm::C3 && get(bondee) == idatm::C1 && heavys_[bondee] == 2) {
        // Not the current atom
        auto bondee_bondee = bondee[0] == atom ?
                             bondee[1] : bondee[0];

        // Change current atom if needed!
         if (get(bondee_bondee) != idatm::C1) {
            return idatm::C1;
        }
    }

    // change 2-substituted Npl bound to isolated C2 to N2
    if (get(atom) == idatm::C2 && get(bondee) == idatm::Npl && heavys_[bondee] == 3) {
        if ((bondee[0] == atom || get(bondee[0]) == idatm::C3) &&
            (bondee[1] == atom || get(bondee[1]) == idatm::C3) &&
            (bondee[2] == atom || get(bondee[2]) == idatm::C3)) {

            return idatm::N2;
        }
    }

    // change 1-substituted Npl (terminal -N=N=N) to N1 (e.g. azide)
    if ((get(atom) == idatm::Npl || get(atom) == idatm::N2)
        && heavys_[bondee] == 2 && mol_[bondee].atomic_number() == 7) {

        auto bondee_bondee = bondee[0] == atom ?
                             bondee[1] : bondee[0];

        if (bondee_bondee.atomic_number() == 7) {
            return idatm::N1;
        }
    }

    // Change sulfonyl
    if (get(atom) == idatm::N3 && get(bondee) == idatm::Son) {
        return idatm::Npl;
    }

    // terminal sulfur on tetrahedral center (thiophosphate)
    if (atom.atomic_number() == 16 &&
        (get(bondee) == idatm::P3p ||
         get(bondee) == idatm::Pox ||
         get(bondee) == idatm::Pac)) {
        return idatm::S3m;
    }

    return get(atom);
}

size_t IDATM::fix_special_(const AtomVertex& atom) {
    // Correct peptide bond
    if (get(atom) == idatm::C2) {

        // Check for a tautomer!!!!
        auto n2_index = mol_.size(), o3_index = mol_.size(); // size = invalid
        for (auto bondee : atom.neighbors()) {
            if (get(bondee) == idatm::O3) o3_index = bondee;
            if (get(bondee) == idatm::N2) n2_index = bondee;
        }

        if (o3_index != mol_.size() && n2_index != mol_.size()){
            set(o3_index, idatm::O2);
            set(n2_index, idatm::Npl);
        }
    }

    // Middle azide: change 2-substituted N1 to N1+ (e.g. azide)
    if (get(atom) == idatm::N1 && atom.degree() == 2) {
        if (atom[0].atomic_number() == Element::N &&
            atom[1].atomic_number() == Element::N) {
            return idatm::N1p;
        }
    }

    // change 3-substituted Npl (eg first nitrogen in -N(-H)=N=N) to N2+)
    if ((get(atom) == idatm::Npl || get(atom) == idatm::N2) && atom.degree() == 3) {
        auto has_n1 = false, has_h = false;
        for (auto bondee : atom.neighbors()) {
            if (get(bondee) == idatm::N1) has_n1 = true;
            if (get(bondee) == idatm::N1p)has_n1 = true;
            if (bondee.atomic_number() == 1) has_h = true;
        }
        if (has_n1 && has_h) return idatm::N2p;
    }

    // (R-)P(=N-)(=N-). We need to change the Ns from Npl to N2
    if (atom.atomic_number() == 15 && atom.degree() == 3) {
        auto invalid = mol_.size();
        auto n_index1 = invalid, n_index2 = invalid; // size = invalid
        size_t carbon_count = 0;
        for (auto bondee : atom.neighbors()) {
            if (bondee.atomic_number() == 6) carbon_count++;
            if (get(bondee) != idatm::Npl && n_index1 != invalid) {
                n_index1 = bondee;
                continue;
            }
            if (get(bondee) != idatm::Npl && n_index1 == invalid) {
                n_index1 = bondee;
                continue;
            }
        }
        if (carbon_count == 1 && n_index1 != invalid && n_index2 != invalid) {
            set(n_index1, idatm::N2);
            set(n_index2, idatm::N2);
        }
    }

    // (R-)(R-)S(=N-)(=N-). We need to change the Ns from Npl to N2
    if (atom.atomic_number() == 16 && atom.degree() == 4) {
        auto invalid = mol_.size();
        auto n_index1 = invalid, n_index2 = invalid; // size = invalid
        size_t carbon_count = 0;
        for (auto bondee : atom.neighbors()) {
            if (mol_[bondee].atomic_number() == 6) carbon_count++;
            if (get(bondee) != idatm::Npl && n_index1 != invalid) {
                n_index1 = bondee;
                continue;
            }
            if (get(bondee) != idatm::Npl && n_index1 == invalid) {
                n_index1 = bondee;
                continue;
            }
        }
        if (carbon_count == 2 && n_index1 != invalid && n_index2 != invalid) {
            set(n_index1, idatm::N2);
            set(n_index2, idatm::N2);
        }
    }

    // Special ligand POB case
    if (get(atom) == idatm::C2 && atom.degree() == 3) {
        auto has_n3 = false, has_pox = false, has_c3 = false;
        for (auto bondee : atom.neighbors()) {
            if (get(bondee) == idatm::Pox) has_pox = true;
            if (get(bondee) == idatm::C3)  has_c3 = true;
            if (get(bondee) == idatm::N3)  has_n3 = true;
        }

        if (has_c3 && has_n3 && has_pox) return idatm::C3;
    }

    return get(atom);
}

bool IDATM::is_aromatic(size_t atom_id) const {
    switch (get(atom_id)) {
    case idatm::Car: return true; break;
    case idatm::Oar: return true; break;
    case idatm::Oarp: return true; break;
    case idatm::Sar: return true; break;
    case idatm::N2: break; // Need special checks...
    case idatm::N2p: break;
    case idatm::Npl: break;
    case idatm::P: break;
    default:
        return false;
    }

    auto aro_search = aromatic_ring_sizes_.find(atom_id);
    return aro_search != aromatic_ring_sizes_.end();
}

bool IDATM::is_planar(size_t atom_id) const {
    switch (get(atom_id)) {
    // Boron and aluminium
    case idatm::B: case idatm::Al:
    // Carbon
    case idatm::C1:  case idatm::C1m: case idatm::C2:   case idatm::Cac:
    case idatm::Car:
    // Oxygen
    case idatm::Oar: case idatm::Oarp: case idatm::O2:  case idatm::O2m:
    case idatm::O1p:
    // Nitrogen
    case idatm::N2:  case idatm::N2p:  case idatm::Npl: case idatm::Ntr:
    case idatm::N1:  case idatm::N1p:  case idatm::Ngp:
    // Sulfur
    case idatm::Sar: case idatm::S2:
        return true;
    case idatm::P: break; // Aromatic phosphorus is weird, but possible
    default:
        return false;
    }

    auto aro_search = aromatic_ring_sizes_.find(atom_id);
    return aro_search != aromatic_ring_sizes_.end();
}

Hybridization IDATM::hybridization(size_t atom_id) const {
    switch (get(atom_id)) {
    case idatm::C:  case idatm::N:    case idatm::O:    case idatm::S:
    case idatm::P:  case idatm::unk:
        return Hybridization::UNKNOWN;
        break;

    case idatm::C1: case idatm::C1m:  case idatm::N1:   case idatm::N1p:
    case idatm::O1p:
        return Hybridization::SP;
        break;

    case idatm::C2:  case idatm::Car: case idatm::Cac:  case idatm::N2:
    case idatm::N2p: case idatm::Ntr: case idatm::Ngp:
    case idatm::O2:  case idatm::O2m: case idatm::Oarp:
    case idatm::S2:
        return Hybridization::SP2;
        break;

    case idatm::C3:  case idatm::N3:  case idatm::N3p: case idatm::Nox:
    case idatm::Npl: case idatm::O3:  case idatm::O3m: case idatm::Oar:
    case idatm::P3p: case idatm::Pox: case idatm::Pac: case idatm::S3:
    case idatm::S3m: case idatm::S3p: case idatm::Son: case idatm::Sac:
    case idatm::Sxd: case idatm::Sar:
        return Hybridization::SP3;
        break;

    default: return Hybridization::FORCED; break;
    }
}

size_t IDATM::add_atom(size_t idx) {
    switch(mol_[idx].atomic_number()) {
    case Element::H:
        push_back(idatm::H);
        break;
    case Element::C:
        push_back(idatm::C3);
        break;
    case Element::N:
        push_back(idatm::N3);
        break;
    case Element::O:
        push_back(idatm::O3);
        break;
    case Element::S:
        push_back(idatm::S3);
        break;
    case Element::P:
        push_back(idatm::P3p);
        break;
    default:
        push_back(idatm_mask.at(Element::Name[mol_[idx].atomic_number()]));
        break;
    }

    mapped_.push_back(false);
    heavys_.push_back(0);

    return get(idx);
}

void IDATM::add_bond(size_t idx1, size_t idx2, Bond::Order bo) {

    if (mol_[idx1].atomic_number() != Element::H) {
        ++heavys_[idx2];
    }

    if (mol_[idx2].atomic_number() != Element::H) {
        ++heavys_[idx1];
    }

    switch (bo) {
    case Bond::SINGLE:
        if (mol_[idx1].atomic_number() == Element::C &&
            mol_[idx2].atomic_number() == Element::H) {
            set(idx2, idatm::HC);
            mapped_[idx2] = true;
        }
        if (mol_[idx2].atomic_number() == Element::C &&
            mol_[idx1].atomic_number() == Element::H) {
            set(idx1, idatm::HC);
            mapped_[idx1] = true;
        }
        break;
    case Bond::DOUBLE:
    case Bond::TRIPLE:
        set(idx1,
            assignBondOrderType(mol_[idx1].atomic_number(),
                                bo,
                                get(idx1))
        );

        set(idx2,
            assignBondOrderType(mol_[idx2].atomic_number(),
                                bo,
                                get(idx2))
        );
        break;
    default:
        break;
    }
}

template<> std::string Spear::atomtype_name_for_id<IDATM>(size_t id) {
    return idatm_unmask[id];
}

template<> size_t Spear::atomtype_id_for_name<IDATM>(const std::string& name) {
    return idatm_mask.at(name);
}

template<> size_t Spear::atomtype_id_count<IDATM>() {
    assert(idatm_mask.size() == idatm::unk + 1);
    assert(idatm_mask.size() == sizeof(idatm_unmask) / sizeof(char*));
    return idatm_mask.size();
}

template<> double Spear::van_der_waals<IDATM>(size_t id) {
    switch (id) {
    case idatm::Ag: return 1.72;
    case idatm::Ar: return 1.88;
    case idatm::As: return 1.85;
    case idatm::Au: return 1.66;
    case idatm::Br: return 1.85;

    case idatm::C:  case idatm::C1:  case idatm::C1m: case idatm::C2:
    case idatm::C3: case idatm::Cac: case idatm::Car: return 1.7;

    case idatm::Cd: return 1.58;
    case idatm::Cl: return 1.75;
    case idatm::Cu: return 1.4;

    case idatm::H: case idatm::HC:
    case idatm::D: case idatm::DC: return 1.2;

    case idatm::F: return 1.47;
    case idatm::Ga: return 1.87;
    case idatm::He: return 1.4; 
    case idatm::Hg: return 1.55;
    case idatm::I: return 1.98;
    case idatm::In: return 1.93;
    case idatm::K: return 2.75;
    case idatm::Kr: return 2.02;
    case idatm::Li: return 1.82;
    case idatm::Mg: return 1.73;

    case idatm::N:   case idatm::N1:  case idatm::N1p: case idatm::N2:
    case idatm::N2p: case idatm::N3:  case idatm::N3p: case idatm::Ngp:
    case idatm::Nox: case idatm::Npl: case idatm::Ntr: return 1.55;

    case idatm::Na: return 2.27;
    case idatm::Ne: return 1.54;
    case idatm::Ni: return 1.63;

    case idatm::O:   case idatm::O1: case idatm::O1p: case idatm::O2:
    case idatm::O2m: case idatm::O3: case idatm::O3m: case idatm::Oar:
    case idatm::Oarp: return 1.52;

    case idatm::Pb: return 2.02;
    case idatm::Pd: return 1.63;
    case idatm::Pt: return 1.72;

    case idatm::P:   case idatm::P3p: case idatm::Pac: case idatm::Pox:
    case idatm::S:   case idatm::S2:  case idatm::S3:  case idatm::S3m:
    case idatm::S3p: case idatm::Sac: case idatm::Sar: case idatm::Son:
    case idatm::Sxd: return 1.8;

    case idatm::Se: return 1.9;
    case idatm::Si: return 2.1;
    case idatm::Sn: return 2.17;
    case idatm::Te: return 2.06;
    case idatm::Tl: return 1.96;
    case idatm::U: return 1.86;
    case idatm::Xe: return 2.16;
    case idatm::Zn: return 1.39;
    default:
        return 2.0;
    }
}
