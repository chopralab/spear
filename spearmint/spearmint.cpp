#include <memory>

#include "spearmint.h"

#include "chemfiles/Trajectory.hpp"

#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Graph_impl.hpp"
#include "spear/Grid.hpp"
#include "spear/ScoringFunction.hpp"

#include "spear/atomtypes/IDATM.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"

std::unique_ptr<Spear::Molecule> receptor;
std::unique_ptr<Spear::Molecule> ligand;
std::unique_ptr<Spear::ScoringFunction> score;
std::unique_ptr<Spear::Grid> gridrec;
std::string error_string = "";

const char* spear_get_error() {
    return error_string.c_str();
}

void set_error(const std::string& error) {
    error_string = "[Spearmint] " + error;
}

size_t spear_receptor_atom_count() {
    if (receptor == nullptr) {
        error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    return receptor->size();
}

size_t spear_receptor_atoms(float* pos) {
    if (receptor == nullptr) {
        error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        size_t current = 0;
        for (auto& rpos : receptor->positions()) {
            pos[current * 3 + 0] = static_cast<float>(rpos[0]);
            pos[current * 3 + 1] = static_cast<float>(rpos[1]);
            pos[current * 3 + 2] = static_cast<float>(rpos[2]);
        }

        return 1;
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating receptor atom arrays: ") + e.what();
        return 0;
    }
}

size_t spear_receptor_atom_details(char* cids, size_t* resi,
                                   char* resn, size_t* elements) {
    if (receptor == nullptr) {
        error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {

        for (auto residue : receptor->topology().residues()) {
            for (auto atom : residue) {
                cids[atom] = residue.get("chain_name")->as_string()[0]; // FIXME
                resi[atom] = *residue.id();
                resn[atom] = residue.name()[0]; // FIXME
                elements[atom] = static_cast<size_t>((*receptor)[atom].atomic_number());
            }
        }

        return 1;
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating receptor atom arrays: ") + e.what();
        return 0;
    }
}

size_t spear_receptor_bond_count() {
    if (receptor == nullptr) {
        error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {
        return receptor->topology().bonds().size();
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating receptor bond arrays: ") + e.what();
        return 0;
    }
}

size_t spear_receptor_bonds(size_t* bonds) {
    if (receptor == nullptr) {
        error_string = std::string("You must run initialize_receptor first");
        return 0;
    }

    try {

        size_t i = 0;
        auto& bos = receptor->topology().bond_orders();
        for (auto& a : receptor->topology().bonds()) {
            bonds[i * 3 + 0] = a[0];
            bonds[i * 3 + 1] = a[1];
            bonds[i * 3 + 2] = static_cast<size_t>(bos[i]);
            ++i;
        }

        return i;
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating receptor bond arrays: ") + e.what();
        return 0;
    }
}

size_t spear_ligand_atom_count() {
    if (ligand == nullptr) {
        error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    return ligand->size();
}

size_t spear_ligand_atoms(float* pos) {
    if (ligand == nullptr) {
        error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        size_t current = 0;
        for (auto& rpos : ligand->positions()) {
            pos[current * 3 + 0] = static_cast<float>(rpos[0]);
            pos[current * 3 + 1] = static_cast<float>(rpos[1]);
            pos[current * 3 + 2] = static_cast<float>(rpos[2]);
        }

        return 1;
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating ligand atom arrays: ") + e.what();
        return 0;
    }
}

size_t spear_ligand_atom_details(size_t* elements) {
    if (ligand == nullptr) {
        error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {

        for (auto atom : *ligand) {
            elements[atom] = static_cast<size_t>(atom.atomic_number());
        }

        return 1;
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating ligand atom arrays: ") + e.what();
        return 0;
    }
}

size_t spear_ligand_bond_count() {
    if (ligand == nullptr) {
        error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        return ligand->topology().bonds().size();
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating ligand bond arrays: ") + e.what();
        return 0;
    }
}

size_t spear_ligand_bonds(size_t* bonds) {
    if (ligand == nullptr) {
        error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {

        size_t i = 0;
        auto& bos = ligand->topology().bond_orders();
        for (auto& a : ligand->topology().bonds()) {
            bonds[i * 3 + 0] = a[0];
            bonds[i * 3 + 1] = a[1];
            bonds[i * 3 + 2] = static_cast<size_t>(bos[i]);
            ++i;
        }

        return i;
    } catch (std::exception& e) {
        error_string =
            std::string("Error creating ligand bond arrays: ") + e.what();
        return 0;
    }
}


size_t ligand_get_neighbors(size_t atom_idx, size_t* neighbors) {
    if (ligand == nullptr) {
        error_string = std::string("You must run initialize_ligand first");
        return 0;
    }

    try {
        size_t current = 0;
        for (auto neighbor : (*ligand)[atom_idx].neighbors()) {
            neighbors[current++] = static_cast<size_t>(neighbor);
        }
    } catch (std::exception& e) {
        error_string = std::string("Error creating atom arrays");
        return 0;
    }

    return 1;
}

size_t spear_initialize_scoring(const char* data_dir) {
    try {
        auto atomtype_name  = receptor->add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
        auto atomtype_name2 = ligand->add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);

        auto ptypes = receptor->get_atomtype(atomtype_name);
        auto ltypes = ligand->get_atomtype(atomtype_name);

        std::unordered_set<size_t> all_types;
        std::copy(ptypes->cbegin(), ptypes->cend(), std::inserter(all_types, all_types.begin()));
        std::copy(ltypes->cbegin(), ltypes->cend(), std::inserter(all_types, all_types.begin()));

        // Remove hydrogen types
        all_types.erase(47);
        all_types.erase(48);

        std::ifstream csd_distrib((std::string(data_dir) + "/csd_distributions.dat").c_str());
        if (!csd_distrib) {
            set_error("Could not open distribution file.");
            return 0;
        }

        Spear::AtomicDistributions atomic_distrib = Spear::read_atomic_distributions<Spear::IDATM>(csd_distrib);

        using Spear::Bernard12;
        auto options = Bernard12::Options(Bernard12::RADIAL | Bernard12::MEAN | Bernard12::REDUCED);
        score = std::make_unique<Bernard12>(options, 6.0, atomic_distrib, atomtype_name, all_types);
    } catch (const std::exception& e) {
        set_error(std::string("Error in loading ligand ") + e.what());
        return 0;
    }
    return 1;
}

float spear_calculate_score() {
    return score->score(*gridrec, *receptor, *ligand);
}

size_t set_positions_ligand(const size_t* atoms, const float* positions,
                            size_t size) {
}

size_t set_positions_receptor(const size_t* atoms, const float* positions,
                              size_t size) {
}

size_t spear_ligand_is_adjacent(size_t atom1, size_t atom2) {
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }
    auto& graph = ligand->graph();
    auto edge = boost::edge(atom1, atom2, graph);
    return edge.second;
}

size_t spear_add_ligand_bond(size_t atom1, size_t atom2) {
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }

    auto av1 = (*ligand)[atom1];
    if (av1.degree() >= av1.expected_bonds()) {
        return 0;
    }

    auto av2 = (*ligand)[atom2];
    if (av2.degree() >= av2.expected_bonds()) {
        return 0;
    }

    ligand->add_bond(atom1, atom2);
    return 1;
}

size_t remove_ligand_bond(size_t atom1, size_t atom2) {
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }

    return 1;
}
