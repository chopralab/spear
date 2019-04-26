#include "spear/forcefields/GAFF_FF.hpp"
#include "spear/atomtypes/GAFF.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Simulation.hpp"
#include "spear/Molecule_impl.hpp"

#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using Spear::GAFF_FF;
using Spear::GAFF;

static void update_chfl_frame(chemfiles::Frame& frame,
                              const std::vector<Eigen::Vector3d>& new_pos) {
    CHECK(frame.size() == new_pos.size());
    auto pos_span = frame.positions();
    for (size_t i = 0; i < frame.size(); ++i) {
        pos_span[i] = {new_pos[i][0], new_pos[i][1], new_pos[i][2]};
    }
}

TEST_CASE("Read gaff.dat file") {
    std::ifstream gaff_dat("share/gaff_2.1.dat");
    GAFF_FF ff(gaff_dat);

    CHECK(ff.num_atom_types() == 83);
    CHECK(ff.num_bonds() == 840);
    CHECK(ff.num_angles() == 4614);
    CHECK(ff.num_torsions() == 980);
    CHECK(ff.num_impropers() == 35);

    SECTION("Tibolone") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto frame = traj.read();
        auto mol = Spear::Molecule(frame);
        mol.add_hydrogens();

        Spear::Simulation sim;

        CHECK_THROWS(sim.add_molecule(mol, ff));

        mol.add_atomtype<GAFF>();
        sim.add_molecule(mol, ff);

        /*chemfiles::Trajectory otraj("tibolone.sdf.gz", 'w');
        chemfiles::Frame start;
        start.resize(mol.size());
        start.set_topology(mol.topology());
        auto pos = sim.positions();

        while (sim.time() <= 10.0) {
            sim.dynamic_steps(100);
            pos = sim.positions();
            update_chfl_frame(start, pos);
            otraj.write(start);
        }*/

        sim.dynamic_steps(100);
    }

    ff.reset();

    SECTION("3qox") {
        auto ltraj = chemfiles::Trajectory("data/3qox_ligand.sdf");
        auto ligand = Spear::Molecule(ltraj.read());
        auto idatm = ligand.add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
        ligand.set_default_atomtype(idatm);
        ligand.add_atomtype<GAFF>();

        Spear::Simulation sim;
        sim.add_molecule(ligand, ff);

        /*chemfiles::Trajectory otraj("3qox.mol2.gz", 'w');

        auto topo = protein.topology();
        chemfiles::Residue ligres("LIG", 9999);
        for (auto liga : ligand.topology()) {
            topo.add_atom(liga);
            ligres.add_atom(topo.size() - 1);
        }
        for (auto bond : ligand.topology().bonds()) {
            topo.add_bond(bond[0], bond[1]);
        }
        topo.add_residue(ligres);

        chemfiles::Frame start;
        start.resize(topo.size());
        start.set_topology(topo);
        auto pos = sim.positions();

        while (sim.time() <= 10.0) {
            sim.dynamic_steps(100);
            pos = sim.positions();
            update_chfl_frame(start, pos);
            otraj.write(start);
        }*/

        sim.dynamic_steps(100);
    }
}
