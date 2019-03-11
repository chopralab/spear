#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

#include "spear/scoringfunctions/Bernard12.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "spear/Grid.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

TEST_CASE("Protein-Ligand Score") {
    auto ptraj = chemfiles::Trajectory("data/3qox_pocket.pdb");
    auto protein = Molecule(ptraj.read());
    auto atomtype_name = protein.add_atomtype<IDATM>(AtomType::GEOMETRY);

    auto ltraj = chemfiles::Trajectory("data/3qox_ligand.sdf");
    auto ligand = Molecule(ltraj.read());
    auto atomtype_name2 = ligand.add_atomtype<IDATM>(AtomType::GEOMETRY);

    CHECK(atomtype_name == atomtype_name2);

    auto ptypes = protein.get_atomtype(atomtype_name);
    auto ltypes = ligand.get_atomtype(atomtype_name);

    std::unordered_set<size_t> all_types;
    std::copy(ptypes->cbegin(), ptypes->cend(), std::inserter(all_types, all_types.begin()));
    std::copy(ltypes->cbegin(), ltypes->cend(), std::inserter(all_types, all_types.begin()));

    // Remove hydrogen types
    all_types.erase(47);
    all_types.erase(48);

    std::ifstream csd_distrib("share/csd_distributions.dat");

    AtomicDistributions atomic_distrib = read_atomic_distributions<IDATM>(csd_distrib);
    auto options = Bernard12::Options(Bernard12::RADIAL | Bernard12::MEAN | Bernard12::REDUCED);
    Bernard12 scoring_func(options, 6.0, atomic_distrib, atomtype_name, all_types);

    std::vector<Eigen::Vector3d> positions;
    for (auto& pos : protein.frame().positions()) {
        positions.push_back({pos[0], pos[1], pos[2]});
    }

    auto grid = Grid(positions);
    CHECK(std::fabs(scoring_func.score(grid, protein, ligand) - -61.8901) < 1e-3);

    auto junk = chemfiles::Trajectory("data/3qox_pocket.pdb");
    auto junk2 = Molecule(junk.read());
    CHECK_THROWS(scoring_func.score(grid, protein, junk2));
}
