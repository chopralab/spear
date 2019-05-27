#include <memory>

#include "spearmint.h"

#include "chemfiles/Trajectory.hpp"

#include "spear/Molecule.hpp"
#include "spear/Grid.hpp"
#include "spear/ScoringFunction.hpp"

#include "spear/atomtypes/IDATM.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"

std::unique_ptr<Spear::Molecule> receptor;
std::unique_ptr<Spear::Molecule> ligand;
std::unique_ptr<Spear::ScoringFunction> score;
std::unique_ptr<Spear::Grid> gridrec;
static thread_local std::string error_string = "";

const char* spear_get_error() {
    auto cstring = error_string.c_str();
    char* copy = new char[error_string.size() + 1];
    strcpy(copy, cstring);
    return copy;
}

void set_error(const std::string& error) {
    error_string = "[Spearmint] " + error;
}

uint64_t get_atom_positions(const Spear::Molecule& mol, float* pos) {
    try {
        size_t current = 0;
        for (auto& rpos : mol.positions()) {
            pos[current * 3 + 0] = static_cast<float>(rpos[0]);
            pos[current * 3 + 1] = static_cast<float>(rpos[1]);
            pos[current * 3 + 2] = static_cast<float>(rpos[2]);
            current++;
        }

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error creating atom arrays: ") + e.what());
        return 0;
    }
}

uint64_t get_bonds(Spear::Molecule& mol, size_t* bonds) {
    try {
        size_t i = 0;
        auto& bos = mol.topology().bond_orders();
        for (auto& a : mol.topology().bonds()) {
            bonds[i * 3 + 0] = a[0];
            bonds[i * 3 + 1] = a[1];
            bonds[i * 3 + 2] = static_cast<size_t>(bos[i]);
            ++i;
        }

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error setting bond arrays: ") + e.what());
        return 0;
    }
}

uint64_t set_positions(Spear::Molecule& mol, const float* positions) {
    auto size = mol.size();
    std::vector<Eigen::Vector3d> posvector;
    posvector.reserve(size / 3);

    for (size_t i = 0; i < size; ++i) {
        posvector.emplace_back(
            positions[i * 3 + 0],
            positions[i * 3 + 1],
            positions[i * 3 + 2]
        );
    }

    mol.set_positions(std::move(posvector));

    return 1;
}

uint64_t spear_initialize_scoring(const char* data_dir) {
    if (ligand == nullptr || receptor == nullptr) {
        set_error("You must initialize ligand and receptor first");
        return 0;
    }

    try {
        auto atomtype_name  = receptor->add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
        auto atomtype_name2 = ligand->add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);

        receptor->atomtype(atomtype_name);
        ligand->atomtype(atomtype_name);

        std::ifstream csd_distrib((std::string(data_dir) + "/csd_distributions.dat").c_str());
        if (!csd_distrib) {
            set_error("Could not open distribution file.");
            return 0;
        }

        Spear::AtomicDistributions atomic_distrib = Spear::read_atomic_distributions<Spear::IDATM>(csd_distrib);

        using Spear::Bernard12;
        auto options = Bernard12::Options(Bernard12::RADIAL | Bernard12::MEAN | Bernard12::COMPLETE);
        score = std::make_unique<Bernard12>(options, 15.0, atomic_distrib, atomtype_name);
    } catch (const std::exception& e) {
        set_error(std::string("Error in loading ligand ") + e.what());
        return 0;
    }
    return 1;
}

float spear_calculate_score() {
    if (ligand == nullptr || receptor == nullptr) {
        set_error("You must initialize ligand and receptor first");
        return 0.000;
    }

    if (score == nullptr) {
        set_error("You must run initialize_score first.");
        return 0.000;
    }

    return static_cast<float>(score->score(*gridrec, *receptor, *ligand));
}
