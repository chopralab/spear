// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include "spear/forcefields/GAFF_FF.hpp"
#include "spear/atomtypes/GAFF.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Geometry.hpp"

#include "OpenMM.h"

#include <string>
#include <iostream>
#include <regex>

using namespace Spear;

GAFF_FF::GAFF_FF(std::istream& input) {
    try {
        read_dat_file_(input);
    } catch (const std::regex_error& e) {
        std::cerr << e.what() << std::endl;
    }
}

template<typename T>
T& get_force(OpenMM::System& system, int& force_id) {
    T* force = nullptr;
    if (force_id == -1) {
        force = new T();
        force_id = system.addForce(force);
    } else {
        force = dynamic_cast<T*>(&system.getForce(force_id));
    }

    return *force;
}

void GAFF_FF::add_forces(const std::vector<std::reference_wrapper<const Molecule>>& mols, OpenMM::System& system) const {
    
    auto& nonbond = get_force<OpenMM::NonbondedForce>(system, non_bond_force_);

    size_t added = 0;
    std::vector<std::pair<int,int>> bonds;
    for (auto& mol : mols) {
        auto gaff_types = mol.get().atomtype("gaff");
        auto& all_types = gaff_types->as_vec();

        for (auto av : mol.get()) {
            auto lookup = type_to_atom_.find(all_types[av]);
            if (lookup == type_to_atom_.end()) {
                std::cerr << "Warning: atom '"
                          << atomtype_name_for_id<GAFF>(all_types[av]) << "'"
                          << " not found.\n";
            }
            nonbond.addParticle(0.0, lookup->second.sigma, lookup->second.epsilon); //FIXME
        }

        auto& topo = mol.get().topology();
        for (auto bond : topo.bonds()) {
            bonds.emplace_back(static_cast<int>(bond[0] + added),
                               static_cast<int>(bond[1] + added)
            );
        }
        added += mol.get().size();
    }

    nonbond.createExceptionsFromBonds(bonds, 0.5, 0.833333333);
}

void GAFF_FF::add_forces(const Molecule& mol, OpenMM::System& system) const {
    auto gaff_types = mol.atomtype("gaff");

    if (static_cast<int>(mol.size()) > system.getNumParticles()) {
        throw std::runtime_error("You must add all particles to the system first.");
    }
    auto added = static_cast<size_t>(system.getNumParticles()) - mol.size();

    if (gaff_types == nullptr) {
        throw std::runtime_error("You must add the 'GAFF' atomtype to the molecule.");
    }

    auto& all_types = gaff_types->as_vec();

    auto& hbond = get_force<OpenMM::HarmonicBondForce>(system, bond_force_);

    for (auto bond : mol.topology().bonds()) {
        auto lookup = type_to_bond_.find({all_types[bond[0]],
                                          all_types[bond[1]]
        });

        if (lookup == type_to_bond_.end()) {
            std::cerr << all_types[bond[0]] << " " << all_types[bond[1]] << std::endl;
            std::cerr << "Warning: bond '"
                      << atomtype_name_for_id<GAFF>(all_types[bond[0]]) << "-"
                      << atomtype_name_for_id<GAFF>(all_types[bond[1]]) << "'"
                      << " not found.\n";
            continue;
        }

        hbond.addBond(static_cast<int>(bond[0] + added),
                      static_cast<int>(bond[1] + added),
            lookup->second.length, lookup->second.k
        );
    }

    auto& hangle = get_force<OpenMM::HarmonicAngleForce>(system, angle_force_);

    for (auto angle : mol.topology().angles()) {
        auto lookup = type_to_angle_.find({all_types[angle[0]],
                                           all_types[angle[1]],
                                           all_types[angle[2]]
        });

        if (lookup == type_to_angle_.end()) {
            std::cerr << "Warning: angle '"
                      << atomtype_name_for_id<GAFF>(all_types[angle[0]]) << "-"
                      << atomtype_name_for_id<GAFF>(all_types[angle[1]]) << "-"
                      << atomtype_name_for_id<GAFF>(all_types[angle[2]]) << "'"
                      << " not found.\n";
            continue;
        }

        hangle.addAngle(static_cast<int>(angle[0] + added),
                        static_cast<int>(angle[1] + added),
                        static_cast<int>(angle[2] + added),
            lookup->second.theta, lookup->second.k
        );
    }

    auto& torsion = get_force<OpenMM::PeriodicTorsionForce>(system, torsion_force_);

    for (auto dihedral : mol.topology().dihedrals()) {
        auto lookup = type_to_torsion_.find({all_types[dihedral[0]],
                                             all_types[dihedral[1]],
                                             all_types[dihedral[2]],
                                             all_types[dihedral[3]]
        });

        if (lookup == type_to_torsion_.end()) {
            lookup = type_to_torsion_.find({0,
                                            all_types[dihedral[1]],
                                            all_types[dihedral[2]],
                                            0
            });
        }

        if (lookup == type_to_torsion_.end()) {
            std::cerr << "Warning: torison '"
                      << atomtype_name_for_id<GAFF>(all_types[dihedral[0]]) << "-"
                      << atomtype_name_for_id<GAFF>(all_types[dihedral[1]]) << "-"
                      << atomtype_name_for_id<GAFF>(all_types[dihedral[2]]) << "-"
                      << atomtype_name_for_id<GAFF>(all_types[dihedral[3]]) << "'"
                      << " not found.\n";
            continue;
        }

        torsion.addTorsion(static_cast<int>(dihedral[0] + added),
                           static_cast<int>(dihedral[1] + added),
                           static_cast<int>(dihedral[2] + added),
                           static_cast<int>(dihedral[3] + added),
                           static_cast<int>(std::abs(lookup->second.periodicity)),
                           lookup->second.phase,
                           lookup->second.k
        );
    }

    for (auto improper : mol.topology().impropers()) {
        auto lookup = type_to_improper_.find({all_types[improper[0]],
                                             all_types[improper[1]],
                                             all_types[improper[2]],
                                             all_types[improper[3]]
        });

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0,
                                            all_types[improper[1]],
                                            0,
                                            all_types[improper[2]],
            });
        }

        // No warning, these are rare
        if (lookup == type_to_improper_.end()) {
            continue;
        }

        // Note the switch!
        torsion.addTorsion(static_cast<int>(improper[1] + added),
                           static_cast<int>(improper[0] + added),
                           static_cast<int>(improper[2] + added),
                           static_cast<int>(improper[3] + added),
                           static_cast<int>(std::abs(lookup->second.periodicity)),
                           lookup->second.phase,
                           lookup->second.k
        );
    }
}

std::vector<double> GAFF_FF::masses(const Molecule& mol) const {
    auto gaff_types = mol.atomtype("gaff");

    if (gaff_types == nullptr) {
        throw std::runtime_error("You must add the 'GAFF' atomtype to the molecule.");
    }

    std::vector<double> result;
    result.reserve(mol.size());

    for (auto at : *gaff_types) {
        result.push_back(type_to_atom_.at(at).mass);
    }

    return result;
}

void GAFF_FF::read_dat_file_(std::istream& input) {
    std::smatch m;
    for (std::string line; std::getline(input, line); ) {
        if (line.compare(0, 5, "AMBER") == 0) {
            read_atom_types_(input);
            std::getline(input, line); // Junk line
            read_bonds_(input);
            read_angles_(input);
            read_torsions_(input);
            read_impropers_(input);
            continue;
        }
        if (line.empty()) {
            continue;
        }
        if (line.compare(0, 4, "MOD4") == 0 ||
            line.compare(0, 6, "NONBON") == 0) {
            read_lj_(input);
            continue;
        }
        if (line.compare(0, 3, "END") == 0) break;
    }
}

void GAFF_FF::read_atom_types_(std::istream& input) {
    std::smatch m;
    for (std::string line; std::getline(input, line); ) {
        if (line.empty()) {
            return;
        }

        if (!std::regex_search(line, m,
                              std::regex(R"(^(\S{1,2})\s+(\S+)\s+(\S+))"))) {
            //FIXME: Print Warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF atom type." << std::endl;
            continue;
        }

        if (!m[1].matched || !m[2].matched || !m[3].matched) {
            //FIXME: Print Warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF atom type." << std::endl;
            continue;
        }

        auto type = atomtype_id_for_name<GAFF>(m[1].str());

        if (type == 0) {
            std::cerr << "Warning: Spear does not define: " << m[1].str() << std::endl;
        }

        auto mass = std::stod(m[2].str());
        auto polar= std::stod(m[3].str());
        type_to_atom_.insert({type, {mass, polar, 0.0, 0.0} });
    }
}

void GAFF_FF::read_bonds_(std::istream& input) {
    std::smatch m;
    auto r = std::regex(R"(^(\S{1,2})\s*-(\S{1,2})\s+(\S+)\s+(\S+))", std::regex::optimize);
    for (std::string line; std::getline(input, line); ) {
        if (line.empty()) {
            return;
        }

        if (!std::regex_search(line, m, r)) {
            //FIXME Print Warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF bond type." << std::endl;
            continue;
        }


        if (!m[1].matched || !m[2].matched || !m[3].matched || !m[4].matched) {
            //FIXEME Print Warning
            continue;
        }

        auto t1 = atomtype_id_for_name<GAFF>(m[1].str());
        auto t2 = atomtype_id_for_name<GAFF>(m[2].str());
        auto k = std::stod(m[3].str()) * 2 * OpenMM::KJPerKcal *
                           OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm;
        auto length = std::stod(m[4].str()) * OpenMM::NmPerAngstrom;
        type_to_bond_.insert({ {t1, t2}, {length, k} });
    }
}

void GAFF_FF::read_angles_(std::istream& input) {
    std::smatch m;
    auto r = std::regex(R"(^(\S{1,2})\s*-(\S{1,2})\s*-(\S{1,2})\s+(\S+)\s+(\S+))", std::regex::optimize);
    for (std::string line; std::getline(input, line); ) {
        if (line.empty()) {
            return;
        }

        if (!std::regex_search(line, m, r)) {
            std::cerr << "Warning: '" << line << "' does not define a GAFF angle type." << std::endl;
            continue;
        }

        if (!m[1].matched || !m[2].matched || !m[3].matched || !m[4].matched ||
            !m[5].matched) {
            // Print warning
            continue;
        }

        auto t1 = atomtype_id_for_name<GAFF>(m[1].str());
        auto t2 = atomtype_id_for_name<GAFF>(m[2].str());
        auto t3 = atomtype_id_for_name<GAFF>(m[3].str());
        auto k = std::stod(m[4].str()) * 2 * OpenMM::KJPerKcal;
        auto angle = std::stod(m[5].str()) * OpenMM::RadiansPerDegree;
        type_to_angle_.insert({ {t1, t2, t3}, {angle, k} });
    }
}

void GAFF_FF::read_torsions_(std::istream& input) {
    std::smatch m;
    auto r = std::regex(R"(^(\S{1,2})\s*-(\S{1,2})\s*-(\S{1,2})\s*-(\S{1,2})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+))", std::regex::optimize);
    for (std::string line; std::getline(input, line); ) {
        if (line.empty()) {
            return;
        }

        if (!std::regex_search(line, m, r)) {
            //FIXME Print Warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF torsion type." << std::endl;
            continue;
        }

        if (!m[1].matched || !m[2].matched || !m[3].matched || !m[4].matched ||
            !m[5].matched || !m[6].matched || !m[7].matched || !m[8].matched) {
            //FIXME Print warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF torsion type." << std::endl;
            continue;
        }

        auto t1 = atomtype_id_for_name<GAFF>(m[1].str());
        auto t2 = atomtype_id_for_name<GAFF>(m[2].str());
        auto t3 = atomtype_id_for_name<GAFF>(m[3].str());
        auto t4 = atomtype_id_for_name<GAFF>(m[4].str());

        auto idivf = std::stod(m[5].str());
        auto k = (1 / idivf) * std::stod(m[6].str()) * OpenMM::KJPerKcal;
        auto phase = std::stod(m[7].str()) * OpenMM::RadiansPerDegree;
        auto periodicity = std::abs(std::stod(m[8].str()));

        type_to_torsion_.insert({ {t1, t2, t3, t4}, {phase, k, periodicity} });
    }
}

void GAFF_FF::read_impropers_(std::istream& input) {
    std::smatch m;
    auto r = std::regex(R"(^(\S{1,2})\s*-(\S{1,2})\s*-(\S{1,2})\s*-(\S{1,2})\s+(\S+)\s+(\S+)\s+(\S+))", std::regex::optimize);
    for (std::string line; std::getline(input, line); ) {
        if (line.empty()) {
            return;
        }

        if (!std::regex_search(line, m, r)) {
            //FIXME Print Warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF improper type." << std::endl;
            continue;
        }

        if (!m[1].matched || !m[2].matched || !m[3].matched || !m[4].matched ||
            !m[5].matched || !m[6].matched || !m[7].matched) {
            //FIXME Print warning
            continue;
        }

        auto t1 = atomtype_id_for_name<GAFF>(m[1].str());
        auto t2 = atomtype_id_for_name<GAFF>(m[2].str());
        auto t3 = atomtype_id_for_name<GAFF>(m[3].str());
        auto t4 = atomtype_id_for_name<GAFF>(m[4].str());

        auto k = std::stod(m[5].str()) * OpenMM::KJPerKcal;
        auto phase = std::stod(m[6].str()) * OpenMM::RadiansPerDegree;
        auto periodicity = std::abs(std::stod(m[7].str()));

        // Note the swap, GAFF is i-j-C-k, we use i-C-j-k
        type_to_improper_.insert({ {t1, t3, t2, t4}, {phase, k, periodicity} });
    }
}

void GAFF_FF::read_lj_(std::istream& input) {
    std::smatch m;
    auto r = std::regex(R"(^\s{2}(\S{1,2})\s+(\S+)\s+(\S+))", std::regex::optimize);
    for (std::string line; std::getline(input, line); ) {
        if (line.empty()) {
            return;
        }

        if (!std::regex_search(line, m, r)) {
            //FIXME: Print Warning
            std::cerr << "Warning: '" << line << "' does not define a GAFF L-J type." << std::endl;
            continue;
        }

        if (!m[1].matched || !m[2].matched || !m[3].matched) {
            //FIXME: Print Warning
            continue;
        }

        auto type = atomtype_id_for_name<GAFF>(m[1].str());
        auto sigma = std::stod(m[2].str()) * OpenMM::NmPerAngstrom *
                                             OpenMM::SigmaPerVdwRadius;
        auto epsilon = std::stod(m[3].str()) * OpenMM::KJPerKcal;

        auto type_iter = type_to_atom_.find(type);
        if(type_iter == type_to_atom_.end()) {
            //FIXME Print Warning
            continue;
        }

        type_iter->second.sigma = sigma;
        type_iter->second.epsilon = epsilon;
    }
}
