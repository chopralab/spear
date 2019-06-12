#include "spear/Molecule.hpp"

#include "spear/forcefields/KnowledgeBasedSF.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/atomtypes/GAFF.hpp"
#include "spear/forcefields/AMBER.hpp"
#include "spear/forcefields/GAFF_FF.hpp"
#include "spear/Simulation.hpp"
#include "spear/Grid.hpp"

#include "openmm/internal/SplineFitter.h"

#include "chemfiles.hpp"
#include "OpenMM.h"

#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

static void update_chfl_frame(chemfiles::Frame& frame,
                              const std::vector<Eigen::Vector3d>& new_pos) {
    CHECK(frame.size() == new_pos.size());
    auto pos_span = frame.positions();
    for (size_t i = 0; i < frame.size(); ++i) {
        pos_span[i] = {new_pos[i][0], new_pos[i][1], new_pos[i][2]};
    }
}

auto IDATM_name = atomtype_name_for_id<IDATM>;

TEST_CASE("Protein-Ligand Score") {
    auto ptraj = chemfiles::Trajectory("data/input.pdb").read();
    ptraj.guess_bonds();
    auto protein = Molecule(ptraj);
    auto at_name = protein.add_atomtype<IDATM>(AtomType::GEOMETRY);

    auto ltraj = chemfiles::Trajectory("data/3qox_ligand.sdf");
    auto ligand = Molecule(ltraj.read());
    ligand.remove_hydrogens();
    auto type = ligand.add_atomtype<IDATM>(AtomType::GEOMETRY);
    ligand.set_default_atomtype(type);
    ligand.add_atomtype<GAFF>();

    auto ptypes = protein.atomtype(at_name);
    auto ltypes = ligand.atomtype(at_name);

    std::unordered_set<size_t> all_types;
    std::copy(ptypes->cbegin(), ptypes->cend(), std::inserter(all_types, all_types.begin()));
    //std::copy(ltypes->cbegin(), ltypes->cend(), std::inserter(all_types, all_types.begin()));

    // Remove hydrogen types
    // all_types.erase(47);
    // all_types.erase(48);

    std::ifstream csd_distrib("share/csd_distributions.dat");
    AtomicDistributions atomic_distrib = read_atomic_distributions<IDATM>(csd_distrib);
    auto options = Bernard12::Options(Bernard12::RADIAL | Bernard12::MEAN | Bernard12::COMPLETE);
    Bernard12 scoring_func(options, 15.0, atomic_distrib, at_name, all_types);

    KnowledgeBasedSF ff(scoring_func, at_name);

    for (auto al1 : all_types) {
        for (auto al2 : all_types) {
            if (al2 > al2) continue;
            double B, rep_idx;
            auto pot = ff.obtain_scores_for_pair_(al1, al2, B, rep_idx);
            std::cout << IDATM_name(al1) << " " << IDATM_name(al2) << " " << B << " " << rep_idx << " ";
            std::cout << van_der_waals<IDATM>(al1) + van_der_waals<IDATM>(al2) << std::endl;
            for (auto p : pot) {
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }

    std::ifstream amber_dat("share/ff10.xml");
    std::ifstream tip3p_dat("share/tip3p.xml");
    AMBER amber(amber_dat);
    amber.add_xml_file(tip3p_dat);

    std::ifstream gaff_dat("share/gaff_2.1.dat");
    GAFF_FF gaff(gaff_dat);

    Spear::Simulation::initialize_plugins();
    Spear::Simulation sim;
    sim.add_molecule(protein, amber);
    //sim.add_molecule(ligand, gaff);
    sim.add_non_bonded_force(ff);
    sim.add_langevin(300, 92, 0.002);
    sim.initialize_context("CUDA");

    chemfiles::Trajectory otraj("sah.pdb", 'w');
    chemfiles::Frame start;
    start.resize(protein.size());

    auto topo_copy = protein.topology();
/*
    chemfiles::Residue res("SAH", 999);
    size_t added_count = 0;
    for (auto atom : ligand.topology()) {
        topo_copy.add_atom(atom);
        res.add_atom(protein.size() + added_count);
        added_count++;
    }
    topo_copy.add_residue(res);

    for (auto bond : ligand.topology().bonds()) {
        topo_copy.add_bond(protein.size() + bond[0],
                           protein.size() + bond[1]);
    }*/
    start.set_topology(topo_copy);
    auto pos = sim.positions();
    update_chfl_frame(start, pos);
    otraj.write(start);

    sim.minimize(1e-3, 1000);
    sim.randomize_velocities(300);

    pos = sim.positions();
    update_chfl_frame(start, pos);
    otraj.write(start);

    while (sim.time() <= 1000) {
        sim.dynamic_steps(100);
        pos = sim.positions();
        update_chfl_frame(start, pos);
        otraj.write(start);
    }
}
