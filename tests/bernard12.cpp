#include "spear/Molecule.hpp"
#include "spear/scoringfunctions/Bernard12.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

TEST_CASE("Protein-Ligand Score") {
    auto ptraj = chemfiles::Trajectory("data/3qox_pocket.pdb");
    auto protein = Molecule(ptraj.read());

    auto ltraj = chemfiles::Trajectory("data/3qox_ligand.sdf");
    auto ligand = Molecule(ltraj.read());

    auto ptypes = IDATM(protein, AtomType::GEOMETRY);
    auto ltypes = IDATM(ligand, AtomType::GEOMETRY);

    std::unordered_set<size_t> all_types;
    std::copy(ptypes.cbegin(), ptypes.cend(), std::inserter(all_types, all_types.begin()));
    std::copy(ltypes.cbegin(), ltypes.cend(), std::inserter(all_types, all_types.begin()));

    // Remove hydrogen types
    all_types.erase(47);
    all_types.erase(48);

    std::ifstream csd_distrib("share/csd_distributions.dat");

    AtomicDistributions atomic_distrib = read_atomic_distributions<IDATM>(csd_distrib);
    Bernard12 scoring_func(Bernard12::Options(Bernard12::RADIAL | Bernard12::MEAN | Bernard12::REDUCED), 6.0, atomic_distrib, all_types);
    CHECK(std::fabs(scoring_func.score(protein, ptypes.all_types(), ligand, ltypes.all_types()) - -61.8901) < 1e-3);
}
