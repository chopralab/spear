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
    virtual void add_forces(const Spear::Molecule& mol, OpenMM::System& system) const override {
        OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
        system.addForce(nonbond);

        for (size_t i = 0; i < mol.size(); ++i) {
            nonbond->addParticle(0.0, 0.3350, 0.996);
        }
    }

    virtual std::vector<double> masses(const Spear::Molecule& mol) const override {
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
    auto mol = Spear::Molecule(std::move(start.clone()));

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
    virtual void add_forces(const Spear::Molecule& mol, OpenMM::System& system) const override {
        OpenMM::GBSAOBCForce* gbsa = new OpenMM::GBSAOBCForce();
        gbsa->setSolventDielectric(80.);
        gbsa->setSoluteDielectric(2.);
        system.addForce(gbsa);

        OpenMM::NonbondedForce* nonbond = new OpenMM::NonbondedForce();
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

    virtual std::vector<double> masses(const Spear::Molecule& mol) const override {
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
    auto mol = Spear::Molecule(std::move(start.clone()));

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
