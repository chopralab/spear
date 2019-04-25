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

    SECTION("3qox") {
        auto traj = chemfiles::Trajectory("data/smd_ini.pdb");
        auto frame = traj.read();

        auto& first_res = *frame.topology().residues().begin();
        const_cast<chemfiles::Residue&>(first_res).set("is_n_terminal", true);

        auto& last_res = *(frame.topology().residues().end() - 1);
        const_cast<chemfiles::Residue&>(last_res).set("is_c_terminal", true);

        auto pocket = Spear::Molecule(frame);

        std::cout << pocket.topology().residues()[0].get<chemfiles::Property::BOOL>("is_n_terminal").value_or(false) << std::endl;

        Spear::Simulation sim;
        sim.add_molecule(pocket, ff);

        chemfiles::Trajectory otraj("pocker.pdb.gz", 'w');
        chemfiles::Frame start;
        start.resize(pocket.size());
        start.set_topology(pocket.topology());
        update_chfl_frame(start, pocket.positions());
        otraj.write(start);

        sim.minimize(1e-3, 100);
        auto pos = sim.positions();
        update_chfl_frame(start, pos);
        otraj.write(start);

        auto forces = sim.forces();
        for (auto f : forces) {
            std::cout << f[0] << " " << f[1] << " " << f[2] << std::endl;
        }
return;
        while (sim.time() <= 0.10) {
            sim.dynamic_steps(10);
            pos = sim.positions();
            update_chfl_frame(start, pos);
            otraj.write(start);
        }
    }
}
