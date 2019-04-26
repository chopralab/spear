#include "spear/FunctionalGroup.hpp"
#include "spear/atomtypes/GAFF.hpp"
#include "spear/Molecule_impl.hpp"

#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using Spear::GAFF;
using Spear::atomtype_name_for_id;
using Spear::atomtype_id_for_name;

TEST_CASE("GAFF Types") {
    CHECK(atomtype_name_for_id<GAFF>(18) == "c2");
    CHECK(atomtype_id_for_name<GAFF>("c2") == 18);

    CHECK(atomtype_name_for_id<GAFF>(50) == "n2");
    CHECK(atomtype_id_for_name<GAFF>("n2") == 50);
}

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

    CHECK(atomtype_name_for_id<GAFF>(atomtypes[0]) == "c3");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[1]) == "nh");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[2]) == "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[3]) == "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[4]) == "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[5]) == "nd");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[6]) == "na");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[7]) == "c3");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[8]) == "cc");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[9]) == "c3");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[10])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[11])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[12])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[13])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[14])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[15])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[16])== "nb");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[17])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[30])== "nb"); // complete the ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[18])== "nh");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[19])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[20])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[21])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[22])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[23])== "c3");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[24])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[25])== "ca");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[26])== "s6");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[27])== "n3");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[28])== "o");
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[29])== "o");
}
