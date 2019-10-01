#include "spear/Geometry.hpp"
#include "spear/Molecule.hpp"
#include "spear/Clustering.hpp"
#include <utility>
#include <numeric>
#include <chrono>

#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Greedy representative") {
    chemfiles::Trajectory frog("data/Frog_2.sdf.xz");

    std::vector<Spear::Conformation> confs;
    confs.reserve(frog.nsteps());
    while (!frog.done()) {
        Spear::Molecule mol(frog.read());
        confs.push_back(mol.positions());
    }

    Spear::greedy_remove_non_representative(confs, Spear::square_deviation, 4 * confs[0].size());

    REQUIRE(confs.size() == 110);

    for (size_t i = 0; i < confs.size(); ++i) {
        for (size_t j = i + 1; j < confs.size(); ++j) {
            CHECK(Spear::rmsd(confs[i], confs[j]) > 2);
        }
    }
}
