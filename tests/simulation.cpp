#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Simulation.hpp"
#include "spear/Forcefield.hpp"

#include "chemfiles.hpp"

#include "OpenMM.h"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

class ArLJFluid : public Spear::Forcefield {
public:
    void add_forces(const Spear::Molecule& mol, OpenMM::System& system) const override {
        auto nonbond = new OpenMM::NonbondedForce();
        system.addForce(nonbond);

        for (size_t i = 0; i < mol.size(); ++i) {
            nonbond->addParticle(0.0, 0.3350, 0.996);
        }
    }

    std::vector<double> masses(const Spear::Molecule& mol) const override {
        return std::vector<double>(mol.size(), 39.95);
    }
};

static void update_chfl_frame(chemfiles::Frame& frame,
                              const std::vector<Eigen::Vector3d>& new_pos) {
    CHECK(frame.size() == new_pos.size());
    auto pos_span = frame.positions();
    for (size_t i = 0; i < frame.size(); ++i) {
        pos_span[i] = {new_pos[i][0], new_pos[i][1], new_pos[i][2]};
    }
}

TEST_CASE("Argon") {
    chemfiles::Frame start;
    start.add_atom(chemfiles::Atom("Ar"), {0.0, 0.0, 0.0});
    start.add_atom(chemfiles::Atom("Ar"), {5.0, 0.0, 0.0});
    start.add_atom(chemfiles::Atom("Ar"), {10., 0.0, 0.0});
    auto mol = Spear::Molecule(start);

    Spear::Simulation sim;
    ArLJFluid arlj;
    sim.add_molecule(mol, arlj);
    CHECK(std::abs(sim.time()) < 1e-3); // No time has passed yet

    CHECK(std::abs(-0.6611955602 - sim.potential_energy()) < 1e-3);
    CHECK(std::abs(sim.kinetic_energy()) < 1e-3); // Should be zero

    // Make sure nothing moved
    auto pos = sim.positions();
    CHECK(pos[0] == Eigen::Vector3d{0.0, 0.0, 0.0});
    CHECK(pos[1] == Eigen::Vector3d{5.0, 0.0, 0.0});
    CHECK(pos[2] == Eigen::Vector3d{10., 0.0, 0.0});

    // Make sure nothing is moving
    auto vel = sim.velocities();
    CHECK(vel[0] == Eigen::Vector3d{0.0, 0.0, 0.0});
    CHECK(vel[1] == Eigen::Vector3d{0.0, 0.0, 0.0});
    CHECK(vel[2] == Eigen::Vector3d{0.0, 0.0, 0.0});

    // But they should move soon!
    auto force = sim.forces();
    CHECK(std::abs(force[0][0]) > 1e-3);
    CHECK(std::abs(force[1][0]) < 1e-3); // Forces cancel for this one
    CHECK(std::abs(force[0][0]) > 1e-3);

    chemfiles::Trajectory traj("out.xyz.gz", 'w');
    pos = sim.positions();
    update_chfl_frame(start, pos);
    traj.write(start);
    while (sim.time() <= 10.0) {
        sim.dynamic_steps(10);
        pos = sim.positions();
        update_chfl_frame(start, pos);
        traj.write(start);
    }

    sim.minimize(1e-4, 1000);
    CHECK(std::abs(-2.0229996338 - sim.potential_energy()) < 1e-3);
}

class NaClSystem : public Spear::Forcefield {
public:
    void add_forces(const Spear::Molecule& mol, OpenMM::System& system) const override {
        auto gbsa = new OpenMM::GBSAOBCForce();
        gbsa->setSolventDielectric(80.);
        gbsa->setSoluteDielectric(2.);
        system.addForce(gbsa);

        auto nonbond = new OpenMM::NonbondedForce();
        system.addForce(nonbond);

        const double sigma_na = 1.8680 * OpenMM::NmPerAngstrom * OpenMM::SigmaPerVdwRadius;
        const double epsilon_na = 0.00277 * OpenMM::KJPerKcal;
        const double gbsa_na = 1.992 * OpenMM::NmPerAngstrom;

        const double sigma_cl = 2.4700 * OpenMM::NmPerAngstrom * OpenMM::SigmaPerVdwRadius;
        const double epsilon_cl = 0.1000 * OpenMM::KJPerKcal;
        const double gbsa_cl = 1.735 * OpenMM::NmPerAngstrom;

        for (auto av : mol) {
            if (av.atomic_number() == Spear::Element::Na) {
                nonbond->addParticle(1.0, sigma_na, epsilon_na);
                gbsa->addParticle(1.0, gbsa_na, 0.8);
            }

            if (av.atomic_number() == Spear::Element::Cl) {
                nonbond->addParticle(-1.0, sigma_cl, epsilon_cl);
                gbsa->addParticle(-1.0, gbsa_cl, 0.8);
            }
        }
    }

    std::vector<double> masses(const Spear::Molecule& mol) const override {
        std::vector<double> ret;
        ret.reserve(mol.size());

        for (auto av : mol) {
            if (av.atomic_number() == Spear::Element::Na) {
                ret.push_back(22.99);
            }

            if (av.atomic_number() == Spear::Element::Cl) {
                ret.push_back(35.45);
            }
        }

        return ret;
    }
};

TEST_CASE("NaCl") {
    chemfiles::Frame start;
    start.add_atom(chemfiles::Atom("Na"), { 8,  0,   0});
    start.add_atom(chemfiles::Atom("Cl"), {-8,  0,   0});
    start.add_atom(chemfiles::Atom("Na"), { 0,  9,   0});
    start.add_atom(chemfiles::Atom("Cl"), { 0, -9,   0});
    start.add_atom(chemfiles::Atom("Na"), { 0,  0, -10});
    start.add_atom(chemfiles::Atom("Cl"), { 0,  0,  10});
    auto mol = Spear::Molecule(start);

    Spear::Simulation sim;
    NaClSystem nacl;
    sim.add_molecule(mol, nacl);
    sim.add_langevin(300, 91, 0.002);

    chemfiles::Trajectory traj("nacl.xyz.gz", 'w');
    auto pos = sim.positions();

    traj.write(start);

    while (sim.time() <= 100.0) {
        sim.dynamic_steps(50);
        pos = sim.positions();
        update_chfl_frame(start, pos);
        traj.write(start);
    }
}

class EthaneSystem : public Spear::Forcefield {
public:
    void add_forces(const Spear::Molecule& mol,
                    OpenMM::System& system) const override {

        auto nonbond = new OpenMM::NonbondedForce();
        system.addForce(nonbond);

        const double sigma_h   = 1.4870 * OpenMM::NmPerAngstrom * OpenMM::SigmaPerVdwRadius;
        const double epsilon_h = 0.0157 * OpenMM::KJPerKcal;

        const double sigma_c   = 1.9080 * OpenMM::NmPerAngstrom * OpenMM::SigmaPerVdwRadius;
        const double epsilon_c = 0.1094 * OpenMM::KJPerKcal;

        for (auto av : mol) {
            if (av.atomic_number() == Spear::Element::H) {
                nonbond->addParticle(0.0605, sigma_h, epsilon_h);
            }

            if (av.atomic_number() == Spear::Element::C) {
                nonbond->addParticle(-.1815, sigma_c, epsilon_c);
            }
        }

        auto hbond = new OpenMM::HarmonicBondForce();
        system.addForce(hbond);

        std::vector<std::pair<int,int>> bonds;

        for (auto bond : mol.topology().bonds()) {
            bonds.emplace_back(static_cast<size_t>(bond[0]),
                               static_cast<size_t>(bond[1])
            );
            if (mol[bond[0]].atomic_number() == Spear::Element::C &&
                mol[bond[1]].atomic_number() == Spear::Element::C) {
                hbond->addBond(static_cast<int>(bond[0]), static_cast<int>(bond[1]),
                               1.526 * OpenMM::NmPerAngstrom,
                               310.0 * 2 * OpenMM::KJPerKcal *
                               OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm
                );
            }

            if (mol[bond[0]].atomic_number() == Spear::Element::C &&
                mol[bond[1]].atomic_number() == Spear::Element::H) {
                hbond->addBond(static_cast<int>(bond[0]), static_cast<int>(bond[1]),
                               1.09 * OpenMM::NmPerAngstrom,
                               340.0 * 2 * OpenMM::KJPerKcal *
                               OpenMM::AngstromsPerNm * OpenMM::AngstromsPerNm
                );
            }
        }

        nonbond->createExceptionsFromBonds(bonds, 0.5, 0.5);

        auto hangle = new OpenMM::HarmonicAngleForce();
        system.addForce(hangle);

        for (auto angle : mol.topology().angles()) {
            if (mol[angle[0]].atomic_number() == Spear::Element::C &&
                mol[angle[2]].atomic_number() == Spear::Element::H) {
                hangle->addAngle(static_cast<int>(angle[0]),
                                 static_cast<int>(angle[1]),
                                 static_cast<int>(angle[2]),
                                 109.5 * OpenMM::RadiansPerDegree,
                                 50.0 * 2 * OpenMM::KJPerKcal
                );
            }

            if (mol[angle[0]].atomic_number() == Spear::Element::H &&
                mol[angle[2]].atomic_number() == Spear::Element::H) {
                hangle->addAngle(static_cast<int>(angle[0]),
                                 static_cast<int>(angle[1]),
                                 static_cast<int>(angle[2]),
                                 109.5 * OpenMM::RadiansPerDegree,
                                 35.0 * 2 * OpenMM::KJPerKcal
                );
            }
        }

        auto torsion = new OpenMM::PeriodicTorsionForce();
        system.addForce(torsion);

        for (auto dihedral : mol.topology().dihedrals()) {
            torsion->addTorsion(static_cast<int>(dihedral[0]),
                                static_cast<int>(dihedral[1]),
                                static_cast<int>(dihedral[2]),
                                static_cast<int>(dihedral[3]),
                                3,
                                0 * OpenMM::RadiansPerDegree,
                                0.150 * OpenMM::KJPerKcal
            );
        }
    }

    std::vector<double> masses(const Spear::Molecule& mol) const override {
        std::vector<double> ret;
        ret.reserve(mol.size());
        for (auto av : mol) {
            if (av.atomic_number() == Spear::Element::H) {
                ret.push_back( 1.008);
            }

            if (av.atomic_number() == Spear::Element::C) {
                ret.push_back(12.011);
            }
        }
        return ret;
    }
};

TEST_CASE("Ethane") {
    chemfiles::Frame start;
    start.add_atom(chemfiles::Atom("C"), { -.7605,   0,   0});
    start.add_atom(chemfiles::Atom("C"), {  .7605,   0,   0});
    start.add_atom(chemfiles::Atom("H"), {-1.135, 1.03,   0});
    start.add_atom(chemfiles::Atom("H"), {-1.135, -.51, .89});
    start.add_atom(chemfiles::Atom("H"), {-1.135, -.51,-.89});
    start.add_atom(chemfiles::Atom("H"), { 1.135, 1.03,   0});
    start.add_atom(chemfiles::Atom("H"), { 1.135, -.51, .89});
    start.add_atom(chemfiles::Atom("H"), { 1.135, -.51,-.89});
    start.add_bond(0, 1);
    start.add_bond(0, 2);
    start.add_bond(0, 3);
    start.add_bond(0, 4);
    start.add_bond(1, 5);
    start.add_bond(1, 6);
    start.add_bond(1, 7);
    auto mol = Spear::Molecule(start);

    Spear::Simulation sim;
    EthaneSystem eth;
    sim.add_molecule(mol, eth);
    chemfiles::Trajectory traj("ethane.xyz.gz", 'w');
    auto pos = sim.positions();
    traj.write(start);

    while (sim.time() <= 100.0) {
        sim.dynamic_steps(100);
        pos = sim.positions();
        update_chfl_frame(start, pos);
        traj.write(start);
    }
}
