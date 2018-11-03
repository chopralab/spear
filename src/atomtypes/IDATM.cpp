#include "spear/atomtypes/IDATM.hpp"
#include <map>
#include <string>
#include <locale>

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
    {"HOH",{{"O","O3"}}},
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


enum idatm {
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

const std::unordered_map<uint64_t, std::unordered_map<size_t, idatm>> bond_maps = {
    {6, {{chemfiles::Bond::DOUBLE, idatm::C2}, {chemfiles::Bond::TRIPLE, idatm::C1},
         {chemfiles::Bond::AROMATIC, idatm::Car}}},
    {7, {{chemfiles::Bond::DOUBLE, idatm::N2}, {chemfiles::Bond::TRIPLE, idatm::N1},
         {chemfiles::Bond::AROMATIC, idatm::Npl}, {chemfiles::Bond::AMIDE, idatm::Npl}}},
    {8, {{chemfiles::Bond::DOUBLE, idatm::O2}, {chemfiles::Bond::TRIPLE, idatm::O1},
         {chemfiles::Bond::AROMATIC, idatm::Oar}}},
};

static bool check_dihedrals_for_planarity(const Molecule& mol, const std::set<size_t>& ring) {
    auto& dihedrals = mol.frame().topology().dihedrals();
    // test dihedral angles
    for (const auto& dihedral : dihedrals) {

        if (ring.count(dihedral[0]) == 0 || ring.count(dihedral[1]) == 0 ||
            ring.count(dihedral[2]) == 0 || ring.count(dihedral[3]) == 0) {
            continue;
        }

        auto dangle = mol.frame().dihedral(dihedral[0], dihedral[1],
                                           dihedral[2], dihedral[3]);

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
    const std::vector<EdgeDescriptor> ring_bonds = mol.get_bonds_in(ring);
    for (const auto& bond : ring_bonds) {
        const auto atom1 = boost::source(bond, mol.graph());
        const auto atom2 = boost::target(bond, mol.graph());
        const auto e1 = mol[atom1].atomic_number();
        const auto e2 = mol[atom2].atomic_number();
        double d = mol.frame().distance(atom1, atom2);
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
    double homas = 1.0 - (98.89 / bonds) * sum;
    bool is_aromatic = true;
    if (homas < 0.271) {
        is_aromatic = check_dihedrals_for_planarity(mol, ring);
    }
    return is_aromatic;
}

static size_t freeOxygens(const Molecule& mol, const Spear::AtomVertex& atom,
                          const std::vector<size_t>& heavys) {
    size_t freeOxygens = 0;
    for (auto neighbor : atom) {
        if (mol[neighbor].type() == "O" && heavys[neighbor] == 1) {
            ++freeOxygens;
        }
    }

    return freeOxygens;
}

static size_t assignBondOrderType(uint64_t atomic_number, chemfiles::Bond::BondOrder bo,
                                  size_t current_type) {
    if (current_type == idatm::N2 && bo == chemfiles::Bond::DOUBLE) {
        return idatm::N1;
    }

    // Allenes 
    if (current_type == idatm::C2 && bo == chemfiles::Bond::DOUBLE) {
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

std::vector<size_t> IDATM::type_atoms_3d(const Molecule& mol) {
    // initialize idatm type in Atoms
    atom_types_ = std::vector<size_t>(mol.size(), (size_t)idatm::unk);

    auto indices = boost::get(boost::vertex_index, mol.graph());

    // number of heavy atoms bonded
    heavys_ = std::vector<size_t>(mol.size(), 0);

    // Have we typed the atom?
    mapped_ = std::vector<bool>(mol.size(), false);

    infallible_(mol);
    auto redo = valence_(mol);
    terminal_(mol, redo);
    redo_(mol, redo);
    fix_C2_(mol);
    charges_(mol);
    aromatic_(mol);
    pass8_(mol);
    pass9_(mol);
    pass10_(mol);
    pass11_(mol);

    return atom_types_;
}

std::vector<size_t> IDATM::type_atoms_order(const Molecule& mol) {
    const auto& bond_orders = mol.frame().topology().bond_orders();
    const auto& bonds = mol.frame().topology().bonds();

    atom_types_ = std::vector<size_t>(mol.size(), (size_t)idatm::unk);
    heavys_ = std::vector<size_t>(mol.size(), 0);
    mapped_ = std::vector<bool>(mol.size(), false);

    infallible_(mol);

    for (const auto atom : mol) {
        auto index = atom.index();
        if (mapped_[index]) continue;

        auto freeOs = freeOxygens(mol, atom, heavys_);
        auto valence = atom.neighbor_count();

        switch (atom.atomic_number()) {
            case 6:
                atom_types_[index] = idatm::C3;
                break;
            case 7:
                if (valence == 4) {
                    atom_types_[index] = freeOs >= 1 ? idatm::Nox
                                                     : idatm::N3p;
                } else if (valence == 3) {
                    atom_types_[index] = freeOs >= 2 ? idatm::Ntr
                                                     : idatm::N3;
                } else {
                    atom_types_[index] = idatm::N3;
                }
                break;
            case 8:
                atom_types_[index] = idatm::O3;
                break;
            case 15:
                if (valence == 4) {
                    if (freeOs >= 2) { // phostphate
                        atom_types_[index] = idatm::Pac;
                    } else if (freeOs == 1) { // P-oxide
                        atom_types_[index] = idatm::Pox;
                    } else { // Formally positive SP3 phospohrus
                        atom_types_[index] = idatm::P3p;
                    }
                } else {
                    atom_types_[index] = idatm::P;
                }
                break;
            case 16:
                if (valence == 4) {
                    if (freeOs >= 3) { // Sulfate
                        atom_types_[index] = idatm::Sac;
                    } else if (freeOs >= 1) { // Sulfone
                        atom_types_[index] = idatm::Son;
                    } else { // We don't know! Other!
                        atom_types_[index] = idatm::S;
                    }
                } else if (valence == 3) {
                    atom_types_[index] = freeOs > 0 ? idatm::Sxd
                                                    : idatm::S3p;
                } else {
                    atom_types_[index] = idatm::S3;
                }
                break;
            default:
                break;
        }
    }

    size_t aromatic_count = 0;
    for ( size_t i = 0; i < bonds.size(); ++i) {
        auto bo = bond_orders[i];
        if (bo == chemfiles::Bond::AROMATIC) {
            ++aromatic_count;
        }

        auto index1 = bonds[i][0];
        atom_types_[index1] = assignBondOrderType(mol[index1].atomic_number(),
                                                  bo,
                                                  atom_types_[index1]);

        auto index2 = bonds[i][1];
        atom_types_[index2] = assignBondOrderType(mol[index2].atomic_number(),
                                                  bo,
                                                  atom_types_[index2]);
    }

    fix_C2_(mol);
    charges_(mol);

    if (aromatic_count < 5) {
        aromatic_(mol);
    }

    for (auto& bond : bonds) {
        if (atom_types_[bond[0]] == idatm::N3 &&
            (atom_types_[bond[1]] == idatm::Car)) {
            atom_types_[bond[0]] = idatm::Npl;
        }

        if (atom_types_[bond[1]] == idatm::N3 &&
            (atom_types_[bond[0]] == idatm::Car)) {
            atom_types_[bond[1]] = idatm::Npl;
        }
    }

    pass8_(mol);
    pass9_(mol);
    pass10_(mol);
    pass11_(mol);

    return atom_types_;
}

void IDATM::infallible_(const Molecule& mol) {
    std::locale loc;

    for (const auto atom : mol) {
        if (atom.atomic_number() == 1) { // Hydrogen and deuterium
            bool bondedToCarbon = false;
            for (auto neighbor : atom) {
                if (mol[neighbor].atomic_number() == 6) {
                    bondedToCarbon = true;
                    break;
                }
            }

            const auto& type = atom.type();
            bool isHyd = (std::toupper(type[0], loc) == 'H');

            atom_types_[atom.index()] = bondedToCarbon ? (isHyd ? HC: DC)
                                                       : (isHyd ? H : D);

            mapped_[atom.index()] = true;

            continue;
        }

        size_t heavyCount = 0;
        for (auto neighbor : atom) {
            if (mol[neighbor].atomic_number() > 1) heavyCount++;
        }

        heavys_[atom.index()] = heavyCount;
    }

    // Use templates for "infallible" typing of standard residues
    for (auto& residue : mol.frame().topology().residues()) {
        auto it = standard_residues.find(residue.name());
        if (it != standard_residues.end()) {
            for (auto i : residue) {
                // Hydrogens are mapped
                if (!mapped_[i]) {
                    continue;
                }

                mapped_[i] = true;
                auto a = mol[i];

                // is it the N-terminal residue ?
                if (a.name() == "N" && heavys_[i] == 1) {
                    atom_types_[i] = idatm::N3p;
                    continue;
                }

                // is it C-terminal ?
                if (a.name() == "C" && freeOxygens(mol, a, heavys_) == 2) {
                    atom_types_[i] = idatm::Cac;
                    continue;
                }

                auto it2 = it->second.find(a.name());
                if (it2 == it->second.end()) {
                    throw std::logic_error(
                        "die : cannot find atom name " + a.name() +
                        " of template residue " + residue.name()
                    );
                }
                atom_types_[i] = idatm_mask.at(it2->second);
            }
        }
    }

    for (auto atom : mol) {
        auto index = atom.index();

        // Hydrogens, standard residues
        if (mapped_[index]) {
            continue;
        }

        // undifferentiated types
        auto element = atom.atomic_number();
        if ((element >= 2 && element <= 4) ||
            (element >= 10 && element <= 14) ||
            element >= 17) {
            atom_types_[index] = idatm_mask.at(atom.type());
            mapped_[index] = true; // infallible type
            continue;
        }
    }

}

std::vector<size_t> IDATM::valence_(const Molecule& mol) {
    std::vector<size_t> redo(mol.size());
    for (auto atom : mol) {
        auto index = atom.index();

        // Hydrogens, standard residues
        if (mapped_[index]) {
            continue;
        }

        auto element = atom.atomic_number();
        auto valence = atom.neighbor_count();

        auto freeOs = freeOxygens(mol, atom, heavys_);
        if (valence == 4) { // assume tetrahedral
            if (element == 6) {
                atom_types_[index] = idatm::C3; // must be sp3 carbon
            } else if (element == 7) {
                atom_types_[index] = freeOs >= 1 ? idatm::Nox
                                                 : idatm::N3p;
            } else if (element == 15) {
                if (freeOs >= 2) { // phostphate
                    atom_types_[index] = idatm::Pac;
                } else if (freeOs == 1) { // P-oxide
                    atom_types_[index] = idatm::Pox;
                } else { // Formally positive SP3 phospohrus
                    atom_types_[index] = idatm::P3p;
                }
            } else if (element == 16) {
                if (freeOs >= 3) { // Sulfate
                    atom_types_[index] = idatm::Sac;
                } else if (freeOs >= 1) { // Sulfone
                    atom_types_[index] = idatm::Son;
                } else { // We don't know! Other!
                    atom_types_[index] = idatm::S;
                }
            }
        } else if (valence == 3) {
            auto avgAngle = 0.0;
            for (size_t n1 = 0; n1 < 3; ++n1) {
                for (size_t n2 = n1 + 1; n2 < 3; ++n2) {
                    avgAngle += mol.frame().angle(atom[n1], index, atom[n2]);
                }
            }
            avgAngle /= 3.0;
            avgAngle *= 180 / 3.14149;

            if (element == 6) {
                if (avgAngle < angle23val1) { // angle significantly < 120?
                    atom_types_[index] = idatm::C3; // Then tetrahedral
                } else { // Most likely trigonal planar (some expceptions)
                    atom_types_[index] = freeOs >= 2 ? idatm::Cac
                                                     : idatm::C2;
                }
            } else if (element == 7) {
                if (avgAngle < angle23val1) {
                    atom_types_[index] = idatm::N3; // likely tetrahedral
                } else { 
                    atom_types_[index] = freeOs >= 2 ? idatm::Ntr
                                                     : idatm::Npl;
                }
            } else if (element == 16) { // sulfoxide or formally positive sp3 S
                 atom_types_[index] = freeOs > 0 ? idatm::Sxd
                                                 : idatm::S3p;
            }
        } else if (valence == 2) {
            double ang = mol.frame().angle(atom[0], index, atom[1]) * 180 / 3.14159;
            if (element == 6) {
                if (ang < angle23val1) { // could be tetralhedral, let's redo
                    atom_types_[index] = idatm::C3;
                    redo[index] = 1;
                } else if (ang < angle12val) { // Signficantly less than 180o
                    atom_types_[index] = idatm::C2;
                    if (ang < angle23val2) { // Angle is too small for comfort
                        redo[index] = 3;
                    }
                } else { // It's near 180o, so sp carbon
                    atom_types_[index] = idatm::C1;
                }
            } else if (element == 7) {
                if (ang < angle23val1) { // Check if its tetralhedral, redo
                    atom_types_[index] = idatm::N3;
                    redo[index] = 2;
                } else { // Is it near 180o?
                    atom_types_[index] = ang < angle12val ? idatm::Npl
                                                          : idatm::N1;
                }
            } else if (element == 8) { // Valence of two = tetrahedral
                atom_types_[index] = idatm::O3;
            } else if (element == 16) {
                atom_types_[index] = idatm::S3;
            }
        }

        // For valence 1, ensure that a default type gets assigned.
        if (atom_types_[index] == idatm::unk) {
            atom_types_[index] = idatm_mask.at(atom.type());
        }
    }

    return redo;
}

void IDATM::terminal_(const Molecule& mol, std::vector<size_t>& redo) {
    for (auto atom : mol) {
        auto neighbors = atom.neighbors();
        auto neighbor_count = std::distance(neighbors.first, neighbors.second);

        if (neighbor_count != 1) continue; // only terminal atoms!

        auto bondee = mol[atom[0]];
        auto len = mol.frame().distance(atom.index(), bondee.index());
        auto bondeeType = atom_types_[bondee.index()];
        auto index = atom.index();

        if (atom_types_[index] == idatm::C) { // Default Carbon
            if (mapped_[index]) continue;
            if (len <= p3c1c1 && bondeeType == idatm::C1) { // Check length
               atom_types_[index] = idatm::C1;
            } else if (len <= p3c2c && bondee.atomic_number() == 6) {
                atom_types_[index] = idatm::C2;
            } else if (len <= p3c2n && bondee.atomic_number() == 7) {
                atom_types_[index] = idatm::C2;
            } else {
                atom_types_[index] = idatm::C3;
            }
        } else if (atom_types_[index] == idatm::N) {
            if (mapped_[index]) continue;
            if (len <= p3n1c1 && bondeeType == idatm::C1) {
                atom_types_[index] = idatm::N1;
            } else if (len > p3n3c &&
                       (bondeeType == idatm::C2 ||
                        bondeeType == idatm::C3)) {
                atom_types_[index] = idatm::N3;
            } else if ((len > p3n3n3 && bondeeType == idatm::N3) ||
                       (len > p3n3n2 && bondeeType == idatm::Npl)) {
                atom_types_[index] = idatm::N3;
            } else {
                atom_types_[index] = idatm::Npl;
            }
        } else if (atom_types_[index] == idatm::O) {
            if (bondeeType == idatm::Cac ||
                bondeeType == idatm::Pac ||
                bondeeType == idatm::Sac ||
                bondeeType == idatm::Ntr) {
                if (!mapped_[index]) atom_types_[index] = idatm::O2m;
            } else if (bondeeType == idatm::Nox ||
                       bondeeType == idatm::Pox ||
                       bondeeType == idatm::Son ||
                       bondeeType == idatm::Sxd) {
                if (!mapped_[index]) atom_types_[index] = idatm::O2;
            } else if (len <= p3o2c2 && bondee.atomic_number() == 6) {
                if (!mapped_[index]) atom_types_[index] = idatm::O2;
                if (!mapped_[bondee.index()])
                    atom_types_[bondee.index()] = idatm::C2;
                redo[bondee.index()] = 0; // Don't redo bondee, we are sure.
            } else if (len <= p3o2as && bondee.type() == "As") {
                if (!mapped_[index]) atom_types_[index] = idatm::O2;
            } else {
                if (!mapped_[index]) atom_types_[index] = idatm::O3;
            }
        } else if (atom_types_[index] == idatm::S) {
            if (bondee.atomic_number() == 15) {
                if (!mapped_[index]) atom_types_[index] = idatm::S2;
            } else if (len <= p3s2c2 && bondee.atomic_number() == 6) {
                if (!mapped_[index]) atom_types_[index] = idatm::S2;
                if (!mapped_[bondee.index()])
                    atom_types_[bondee.index()] = idatm::C2;
                redo[bondee.index()] = 0; // Don't redo bondee, we are sure.
            } else if (len <= p3s2as && bondee.type() == "As") {
                if (!mapped_[index]) atom_types_[index] = idatm::S2;
            } else {
                if (!mapped_[index]) atom_types_[index] = idatm::S3;
            }
        }
    }
}

void IDATM::redo_(const Molecule& mol, const std::vector<size_t>& redo) {
    for (auto atom : mol) {
        size_t index = atom.index();
        if (mapped_[index]) continue;
        if (redo[index] == 0) continue;

        bool c3able = false;
        for (auto bondee : atom) {
            auto len = mol.frame().distance(index, bondee);
            auto bondeeElement = *(mol.frame()[bondee].atomic_number());

            if (redo[index] == 1) { // Tetrahedral or planar carbon?
                if ((len > p4c3c && bondeeElement == 6) ||
                    (len > p4c3n && bondeeElement == 7) ||
                    (len > p4c3o && bondeeElement == 8)) {
                    atom_types_[index] = idatm::C3;
                    break; // Done for this redo
                }
                if ((len <= p4c2c && bondeeElement == 6) ||
                    (len <= p4c2n && bondeeElement == 7)) {
                    atom_types_[index] = idatm::C2;
                }
            } else if (redo[index] == 2) { // Tetrahedral or planar nitrogen?
                if ((len <= p4n2c && bondeeElement == 6) ||
                    (len <= p4n2n && bondeeElement == 7)) {
                    atom_types_[index] = idatm::Npl;
                    break; // Done for this redo
                }
            } else {
                if ((len <= p4c2c && bondeeElement == 6) ||
                    (len <= p4c2n && bondeeElement == 7)) {
                    atom_types_[index] = idatm::C2;
                    c3able = false;
                    break;
                }
                if ((len > p4c3c && bondeeElement == 6) ||
                    (len > p4c3n && bondeeElement == 7) ||
                    (len > p4c3o && bondeeElement == 8)) {
                    c3able = true;
                }

                if (len > p4ccnd && bondeeElement == 6) c3able = true;
            }
        }
        if (c3able) atom_types_[index] = idatm::C3;
    }
}

void IDATM::fix_C2_(const Molecule& mol) {
    for (auto atom : mol) {
        auto index = atom.index();
        if (atom_types_[index] != idatm::C2) continue;

        if (freeOxygens(mol, atom, heavys_) >= 2) {
            atom_types_[index] = idatm::Cac;
            continue;
        }

        if (mapped_[index]) continue;

        bool c2possible = false;
        for (auto bondee : atom) {
            if (c2impossibleTypes.count(atom_types_[bondee]) == 0) {
                c2possible = true;
                break;
            }
        }
        if (!c2possible) atom_types_[index] = idatm::C3;
    }
}

void IDATM::charges_(const Molecule& mol) {
    for (auto atom : mol) {
        auto index = atom.index();
        if (atom_types_[index] == idatm::N3) {
            if (mapped_[index]) continue;
            bool positive = true;
            auto c2_index = mol.size();
            for (auto bondee : atom) {
                if (atom_types_[bondee] != idatm::C3 &&
                    atom_types_[bondee] != idatm::H &&
                    atom_types_[bondee] != idatm::D) {
                    positive = false;
                }
                if (atom_types_[bondee] == idatm::C2) {
                    c2_index = bondee;
                }
            }
            if (positive) atom_types_[index] = idatm::N3p;

            if (c2_index != mol.size() && freeOxygens(mol, mol[c2_index], heavys_) == 1) {
                atom_types_[index] = idatm::Npl;
            }
        } else if (atom_types_[index] == idatm::C2) {
            int numNpls = 0;
            bool hasN2 = false;
            for (auto bondee : atom) {
                if ((atom_types_[bondee] == idatm::Npl && !mapped_[bondee]) ||
                    (atom_types_[bondee] == idatm::N3  && heavys_[bondee] == 1 && !mapped_[bondee]) ||
                    atom_types_[bondee] == idatm::Ngp) {
                    // Ng+ possible through template typing
                    numNpls++;
                }
                if (atom_types_[bondee] == idatm::N2 && heavys_[bondee] <= 2 && !mapped_[bondee]) {
                    hasN2 = true;
                }
            }

            bool noplus = false;
            if (numNpls == 3 || (numNpls == 2 && hasN2)) {
                for (auto bondee : atom) {
                    //if (atom_types_[bondee] != idatm::Npl) continue;
                    if (mapped_[bondee]) continue;

                    atom_types_[bondee] = idatm::Ngp;
                    for (auto bondee2 : mol[bondee]) {
                        if ((atom_types_[bondee2] == idatm::C2 ||
                            atom_types_[bondee2] == idatm::Npl) &&
                            bondee2 != index) {
                            atom_types_[bondee] = idatm::Npl;
                            noplus = true;
                            break;
                        }
                    }
                }
            }
            // Reset!
            if (noplus) {
                for (auto bondee : atom) {
                    if (mapped_[bondee]) continue;
                    if (atom_types_[bondee] == idatm::Ngp) {
                        atom_types_[bondee] =  idatm::Npl;
                    }
                }
            }
        } else if (atom_types_[index] == idatm::Cac) {
            for (auto bondee : atom) {
                if (mapped_[bondee]) continue;
                if (mol[bondee].atomic_number() == 8 && heavys_[bondee] == 1) {
                    atom_types_[bondee] =  idatm::O2m;
                }
            }
        }
    }
}

void IDATM::aromatic_(const Molecule& mol) {
    auto rings = mol.rings();
    for (const auto& ring : rings) {
        bool planarTypes = true;
        size_t c3_count = 0;
        for (auto atom : ring) {
            if (possible_aromatic.count(atom_types_[atom]) == 0) {
                if (atom_types_[atom] == idatm::C3 &&
                    mol[atom].neighbor_count() == 2) {
                    ++c3_count;
                    continue;
                }
                // O3/S3/N3/C3 will be changed to Oar/S2/Npl/Car if
                // ring is found to be aromatic
                planarTypes = false;
                break;
            }

            if ((atom_types_[atom] == idatm::C2 ||
                 atom_types_[atom] == idatm::C3 ||
                 atom_types_[atom] == idatm::O3 ||
                 atom_types_[atom] == idatm::S3 ||
                 atom_types_[atom] == idatm::N3) &&
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
                if (atom_types_[atom] == idatm::C2) has_c2 = true;
                if (atom_types_[atom] == idatm::Npl) has_npl = true;
            }
            if (has_c2 && has_npl) continue;
        }

        if (!aromatic(mol, ring)) {
            continue;
        }

        // Correct types to be aromatic
        for (auto atom : ring) {
            aromatic_ring_sizes_.insert({atom, ring.size()});
            if (atom_types_[atom] == idatm::C2) atom_types_[atom] = idatm::Car;
            if (atom_types_[atom] == idatm::C3) atom_types_[atom] = idatm::Car;
            if (atom_types_[atom] == idatm::O3) atom_types_[atom] = idatm::Oar;
            if (atom_types_[atom] == idatm::O2) atom_types_[atom] = idatm::Oarp;
            if (atom_types_[atom] == idatm::S3) atom_types_[atom] = idatm::S2;
            if (atom_types_[atom] == idatm::S2) atom_types_[atom] = idatm::Sar;
            if (atom_types_[atom] == idatm::N3) atom_types_[atom] = idatm::Npl;
        }
    }
}

void IDATM::pass8_(const Molecule& mol) {
    for (auto atom : mol) {
        auto index = atom.index();
        if (mapped_[index]) continue;

        if (atom_types_[index] != idatm::Npl) continue;

        if (heavys_[index] != 2) continue;

        if (atom.neighbor_count() > 2) continue;

        // are both bonded heavy atoms sp3?  If so -> Npl
        bool bothSP3 = true;
        bool bothC = true;
        bool bothN = true;
        bool other = false;
        for (auto bondee : atom) {
            if (sp3hybridized.count(atom_types_[bondee]) == 0) {
                bothSP3 = false;
                // don't break out, need to check both
                // bonded heavies' element types
            }

            if (mol[bondee].atomic_number() == 6) {
                bothN = false;
            } else if (mol[bondee].atomic_number() == 7) {
                bothC = false;
            } else {
                other = true;
                bothN = false;
                bothC = false;
            }
        }
        if (bothSP3) continue;

        if (other) {
            // Metal bond most likely. (C(=O)N:---M), keep Npl
            continue;
        }
        auto avgLen = 0.0;
        for (auto bondee : atom) {
            auto len = mol.frame().distance(index, bondee);
            avgLen += len;
        }
        avgLen /= 2.0;
        if (bothC) {
            if (avgLen <= p8cc2n2h) atom_types_[index] = idatm::N2;
        } else if (bothN) {
            if (avgLen <= p8nn2n2h) atom_types_[index] = idatm::N2;
        } else if (avgLen <= p8cn2n2h) {
            // one N, one C
            atom_types_[index] = idatm::N2;
        }
    }
}

void IDATM::pass9_(const Molecule& mol) {
    // Do two passes
    // First, assign Oar+ and Sar
    for (auto pair : aromatic_ring_sizes_) {
        if (mapped_[pair.first]) continue;

        // aromatic sulfur, always assign??
        if (mol[pair.first].atomic_number() == 16) {
            atom_types_[pair.first] = idatm::Sar;
        }

        // Formally positive oxygen
        if (atom_types_[pair.first] == idatm::Oar && pair.second % 2 == 0) {
            atom_types_[pair.first] = idatm::Oarp;
        }
    }

    // Now we need to retype only Npl
    for (auto pair : aromatic_ring_sizes_) {
        if (mapped_[pair.first]) continue;

        if (atom_types_[pair.first] != idatm::Npl &&
            atom_types_[pair.first] != idatm::N2) continue;

        // change 2-substituted Npl in 6-membered aromatic ring to N2
        // and a 3+-substituted Npl in 6-membered aromatic ring to N2+
        if (pair.second == 6) {
            atom_types_[pair.first] = heavys_[pair.first] == 2 ? idatm::N2 :
                                                                 idatm::N2p;
            continue;
        }

        auto bound_to_car = false;
        auto bound_to_npl = false;
        auto bound_to_oar = false;
        auto bound_to_sar = false; // Note: changed from just S?

        for (auto bondee : mol[pair.first]) {
            switch (atom_types_[bondee]) {
                case idatm::Npl: bound_to_npl = true; break;
                case idatm::Car: bound_to_car = true; break;
                case idatm::Oar: bound_to_oar = true; break;
                case idatm::Oarp: bound_to_oar= true; break;
                case idatm::Sar: bound_to_sar = true; break;
                default: break;
            }
        }

        // change aromatic Npl bound to Npl and Car to N2
        if (bound_to_npl && bound_to_car) {
            atom_types_[pair.first] = idatm::N2;
        }

        // change 2/3-substituted Npl bound to S or O in 5-ring to N2
        if (pair.second == 5 && (bound_to_oar || bound_to_sar) && 
            (heavys_[pair.first] == 2 || heavys_[pair.first] == 3)) {
            atom_types_[pair.first] = idatm::N2;
        }

        // This is odd, we change the first Npl in 5-ring
        // Npl(-H)-Car-Npl(-C,-C,-C) to N2
        // Similar to HIS, but this does not match intuitive atome types...
        if (pair.second == 5 && bound_to_car &&
            (mol[pair.first].neighbor_count() > heavys_[pair.first])) {
            atom_types_[pair.first] = idatm::N2;
        }
    }
}

void IDATM::pass10_(const Molecule& mol) {
    for (auto atom : mol) {
        if (mapped_[atom.index()]) continue;
        if (heavys_[atom.index()] != 1) continue;

        auto index = atom.index();
        auto bondee = atom[0];

        // resonance equivalent terminal oxygen on tetrahedral
        // center (phosphate, sulfate, sulfone...)
        if (atom.atomic_number() == 8 &&
            oxygenConjugatedTetra.count(atom_types_[bondee]) != 0) {
            atom_types_[index] = idatm::O3m;
        }

        // change terminal C3 bound to C1, which in turn is
        // bound to any C but C1 (Ligand ID: 1DJ, Atom Name:
        // CAA) to C1
        if (atom_types_[index] == idatm::C3 &&
            atom_types_[bondee] == idatm::C1 &&
            heavys_[bondee] == 2) {

            // Not the current index
            auto bondee_bondee = mol[bondee][0] == index ?
                                 mol[bondee][1] : mol[bondee][0];

            // Change current index if needed!
            if (atom_types_[bondee_bondee] != idatm::C1) {
                atom_types_[atom.index()] = idatm::C1;
            }
        }

        // change 2-substituted Npl bound to isolated C2 to N2
        if (atom_types_[index] == idatm::C2 &&
            atom_types_[bondee] == idatm::Npl &&
            heavys_[bondee] == 3) {

            if ((mol[bondee][0] == index || atom_types_[mol[bondee][0]] == idatm::C3) &&
                (mol[bondee][1] == index || atom_types_[mol[bondee][1]] == idatm::C3) &&
                (mol[bondee][2] == index || atom_types_[mol[bondee][2]] == idatm::C3)) {

                atom_types_[bondee] = idatm::N2;
            }
        }

        // change 1-substituted Npl (terminal -N=N=N) to N1 (e.g. azide)
        if ((atom_types_[index] == idatm::Npl || atom_types_[index] == idatm::N2)
            && heavys_[bondee] == 2 && mol[bondee].atomic_number() == 7) {

            auto bondee_bondee = mol[bondee][0] == index ?
                                 mol[bondee][1] : mol[bondee][0];

            if (mol[bondee_bondee].atomic_number() == 7) {
                atom_types_[index] = idatm::N1;
            }
        }

        if (atom_types_[index] == idatm::N3 &&
            atom_types_[bondee] == idatm::Son) {
            atom_types_[index] = idatm::Npl;
        }

        // terminal sulfur on tetrahedral center (thiophosphate)
        if (atom.atomic_number() == 16 &&
            (atom_types_[bondee] == idatm::P3p ||
             atom_types_[bondee] == idatm::Pox ||
             atom_types_[bondee] == idatm::Pac)) {
            atom_types_[atom.index()] = idatm::S3m;
        }
    }
}

void IDATM::pass11_(const Molecule& mol) {
    for (auto atom : mol) {
        auto index = atom.index();

        if (mapped_[index]) continue;

        // Correct peptide bond
        if (atom_types_[index] == idatm::C2) {

            // Check for a tautomer!!!!
            auto n2_index = mol.size(), o3_index = mol.size(); // size = invalid
            for (auto bondee : atom) {
                if (atom_types_[bondee] == idatm::O3) o3_index = bondee;
                if (atom_types_[bondee] == idatm::N2) n2_index = bondee;
            }

            if (o3_index != mol.size() && n2_index != mol.size()){
                atom_types_[o3_index] = idatm::O2;
                atom_types_[n2_index] = idatm::Npl;
            }
        }

        // Middle azide: change 2-substituted N1 to N1+ (e.g. azide)
        if (atom_types_[index] == idatm::N1 && atom.neighbor_count() == 2) {
            if (mol[atom[0]].atomic_number() == 7 &&
                mol[atom[1]].atomic_number() == 7) {
                atom_types_[index] = idatm::N1p;
            }
        }

        // change 3-substituted Npl (eg first nitrogen in -N(-H)=N=N) to N2+)
        if ((atom_types_[index] == idatm::Npl || atom_types_[index] == idatm::N2)
            && atom.neighbor_count() == 3) {
            auto has_n1 = false, has_h = false;
            for (auto bondee : atom) {
                if (atom_types_[bondee] == idatm::N1) has_n1 = true;
                if (atom_types_[bondee] == idatm::N1p)has_n1 = true;
                if (mol[bondee].atomic_number() == 1) has_h = true;
            }
            if (has_n1 && has_h) atom_types_[index] = idatm::N2p;
        }

        // (R-)P(=N-)(=N-). We need to change the Ns from Npl to N2
        if (atom.atomic_number() == 15 && atom.neighbor_count() == 3) {
            auto invalid = mol.size();
            auto n_index1 = invalid, n_index2 = invalid; // size = invalid
            size_t carbon_count = 0;
            for (auto bondee : atom) {
                if (mol[bondee].atomic_number() == 6) carbon_count++;
                if (atom_types_[bondee] != idatm::Npl && n_index1 != invalid) {
                    n_index1 = bondee;
                    continue;
                }
                if (atom_types_[bondee] != idatm::Npl && n_index1 == invalid) {
                    n_index1 = bondee;
                    continue;
                }
            }
            if (carbon_count == 1 && n_index1 != invalid && n_index2 != invalid) {
                atom_types_[n_index1] = idatm::N2;
                atom_types_[n_index2] = idatm::N2;
            }
        }

        // (R-)(R-)S(=N-)(=N-). We need to change the Ns from Npl to N2
        if (atom.atomic_number() == 16 && atom.neighbor_count() == 4) {
            auto invalid = mol.size();
            auto n_index1 = invalid, n_index2 = invalid; // size = invalid
            size_t carbon_count = 0;
            for (auto bondee : atom) {
                if (mol[bondee].atomic_number() == 6) carbon_count++;
                if (atom_types_[bondee] != idatm::Npl && n_index1 != invalid) {
                    n_index1 = bondee;
                    continue;
                }
                if (atom_types_[bondee] != idatm::Npl && n_index1 == invalid) {
                    n_index1 = bondee;
                    continue;
                }
            }
            if (carbon_count == 2 && n_index1 != invalid && n_index2 != invalid) {
                atom_types_[n_index1] = idatm::N2;
                atom_types_[n_index2] = idatm::N2;
            }
        }

        // Special ligand POB case
        if (atom_types_[index] == idatm::C2 && atom.neighbor_count() == 3) {
            auto has_n3 = false, has_pox = false, has_c3 = false;
            for (auto bondee : atom) {
                if (atom_types_[bondee] == idatm::Pox) has_pox = true;
                if (atom_types_[bondee] == idatm::C3)  has_c3 = true;
                if (atom_types_[bondee] == idatm::N3)  has_n3 = true;
            }

            if (has_c3 && has_n3 && has_pox) atom_types_[index] = idatm::C3;
        }
    }
}

template<> std::string Spear::atomtype_name<IDATM>() {
    return std::string("IDATM");
}

template<> std::string Spear::atomtype_name_for_id<IDATM>(size_t id) {
    return idatm_unmask[id];
}

template<> size_t Spear::atomtype_id_for_name<IDATM>(std::string name) {
    return idatm_mask.at(name);
}
