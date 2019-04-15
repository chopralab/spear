#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Simulation.hpp"
#include "spear/Forcefield.hpp"

#include "chemfiles.hpp"

#include "OpenMM.h"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

class ArLJFluid : public Spear::Forcefield {
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

    chemfiles::Trajectory traj("out.xyz", 'w');
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
