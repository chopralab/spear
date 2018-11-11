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

    IDATM pidatm, lidatm;
    auto ptypes = pidatm.type_atoms_3d(protein);
    auto ltypes = lidatm.type_atoms_3d(ligand);

    std::unordered_set<size_t> all_types;
    std::copy(ptypes.begin(), ptypes.end(), std::inserter(all_types, all_types.begin()));
    std::copy(ltypes.begin(), ltypes.end(), std::inserter(all_types, all_types.begin()));
    all_types.erase(47);
    all_types.erase(48);

    std::ifstream csd_distrib("share/csd_distributions.dat");

    AtomicDistributions atomic_distrib = read_atomic_distributions<IDATM>(csd_distrib);
    Bernard12 scoring_func(Bernard12::Options(Bernard12::RADIAL | Bernard12::MEAN | Bernard12::REDUCED), 6.0, atomic_distrib, all_types);
    CHECK(std::fabs(scoring_func.score(protein, ptypes, ligand, ltypes) - -61.8901) < 1e-3);
}
