#include "spear/FunctionalGroup.hpp"
#include "spear/atomtypes/GAFF.hpp"
#include "spear/Molecule_impl.hpp"

#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using Spear::GAFF;
using Spear::atomtype_name_for_id;

TEST_CASE("Tibolone") {
    auto traj = chemfiles::Trajectory("data/tibolone.sdf");
    auto mol = Spear::Molecule(traj.read());
    mol.add_hydrogens();
    Spear::GAFF atomtypes(mol);
    for (size_t i = 0; i <= 2; ++i) {    
        CHECK(atomtype_name_for_id<GAFF>(atomtypes[i]) == "c3");
    }
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[3]) == "c2");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[4]) == "c2");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[7]) == "c");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[8]) == "o");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[19]) == "oh");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[20]) == "c1");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[21]) == "c1");
}

TEST_CASE("Pazopanib") {
    auto traj = chemfiles::Trajectory("data/pazopanib.sdf");
    auto mol = Spear::Molecule(traj.read());
    mol.add_hydrogens();
    Spear::GAFF atomtypes(mol);

}
