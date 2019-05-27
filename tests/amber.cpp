#include "spear/forcefields/AMBER.hpp"
#include "spear/Simulation.hpp"
#include "spear/Molecule_impl.hpp"

#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using Spear::AMBER;

static void update_chfl_frame(chemfiles::Frame& frame,
                              const std::vector<Eigen::Vector3d>& new_pos) {
    CHECK(frame.size() == new_pos.size());
    auto pos_span = frame.positions();
    for (size_t i = 0; i < frame.size(); ++i) {
        pos_span[i] = {new_pos[i][0], new_pos[i][1], new_pos[i][2]};
    }
}

TEST_CASE("Read old AMBER file") {
    std::ifstream amber_dat("share/amber10.xml");
    std::ifstream tip3p_dat("share/tip3p.xml");
    AMBER ff(amber_dat);
    ff.add_xml_file(tip3p_dat);

    SECTION("Deca-Alanine") {
        auto traj = chemfiles::Trajectory("data/deca_alanine.pdb");
        auto frame = traj.read();

        auto& first_res = *frame.topology().residues().begin();
        const_cast<chemfiles::Residue&>(first_res).set("is_n_terminal", true);

        auto& last_res = *(frame.topology().residues().end() - 1);
        const_cast<chemfiles::Residue&>(last_res).set("is_c_terminal", true);

        auto pocket = Spear::Molecule(frame);
        pocket.add_bond(0,2);

        Spear::Simulation sim;
        sim.add_molecule(pocket, ff);

        chemfiles::Trajectory otraj("deca_alanine.mmtf", 'w');
        chemfiles::Frame start;
        start.resize(pocket.size());
        start.set_topology(pocket.topology());
        update_chfl_frame(start, pocket.positions());

        sim.minimize(1e-3, 10);
        auto pos = sim.positions();
        update_chfl_frame(start, pos);
        otraj.write(start);

        while (sim.time() <= 5) {
            sim.dynamic_steps(1);
            pos = sim.positions();
            update_chfl_frame(start, pos);
            otraj.write(start);
        }
    }
}

TEST_CASE("Read new AMBER file") {
    std::ifstream amber_dat("share/ff10.xml");
    std::ifstream tip3p_dat("share/tip3p.xml");
    AMBER ff(amber_dat);
    ff.add_xml_file(tip3p_dat);

    Spear::Simulation::initialize_plugins();

    SECTION("Deca-Alanine") {
        auto traj = chemfiles::Trajectory("data/input.pdb.gz");
        auto frame = traj.read();

        auto& first_res = *frame.topology().residues().begin();
        const_cast<chemfiles::Residue&>(first_res).set("is_n_terminal", true);

        auto& last_res = *(frame.topology().residues().begin() + 34);
        const_cast<chemfiles::Residue&>(last_res).set("is_c_terminal", true);

        auto pocket = Spear::Molecule(frame);
        //pocket.add_bond(0,2);

        Spear::Simulation sim;
        ff.set_method(Spear::NonBondedForcefield::Ewald);
        ff.set_cutoff(10.0);
        sim.add_molecule(pocket, ff);
        sim.set_periodic_vectors(frame.cell());
        sim.add_non_bonded_force(ff);
        sim.add_langevin(300.0, 92, 0.002);
        sim.initialize_context("CUDA");

        chemfiles::Trajectory otraj("villin.mmtf", 'w');
        chemfiles::Frame start;
        start.resize(pocket.size());
        start.set_topology(pocket.topology());
        update_chfl_frame(start, pocket.positions());

        sim.minimize(1e-3, 10);
        auto pos = sim.positions();
        update_chfl_frame(start, pos);
        otraj.write(start);
        sim.randomize_velocities(300);

        while (sim.time() <= 20) {
            sim.dynamic_steps(100);
            pos = sim.positions();
            update_chfl_frame(start, pos);
            otraj.write(start);
        }
    }
}