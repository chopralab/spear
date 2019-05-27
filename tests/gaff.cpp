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

    CHECK(atomtype_name_for_id<GAFF>(atomtypes[0]) == "c3"); // Methyl group
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[1]) == "nh"); // briding nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[2]) == "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[3]) == "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[4]) == "ca"); // Pure takes precedence over non-pure
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[5]) == "nd"); // Sp2 nitrogen in inpure ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[6]) == "na"); // Sp3 nitrogen delocalized
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[7]) == "c3"); // Methyl group
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[8]) == "cc"); // Impure arocarbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[9]) == "c3"); // Methyl group
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[10])== "ca"); // Pure takes precedence over non-pure
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[11])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[12])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[13])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[14])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[15])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[16])== "nb"); // Pure ato nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[17])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[30])== "nb"); // Pure aro nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[18])== "nh"); // Briding nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[19])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[20])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[21])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[22])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[23])== "c3"); // methyl group
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[24])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[25])== "ca"); // Aro carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[26])== "s6"); // sulfone sulfur
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[27])== "n3"); // sulfone nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[28])== "o");  // sulfonyl oxygen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[29])== "o");  // sulfonyl oxygen
}

TEST_CASE("SAH (3qox)") {
    auto traj = chemfiles::Trajectory("data/3qox_ligand.sdf");
    auto mol = Spear::Molecule(traj.read());
    Spear::GAFF atomtypes(mol);

    CHECK(atomtype_name_for_id<GAFF>(atomtypes[0]) == "n4"); // Amine
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[1]) == "c3"); // alpha carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[2]) == "c3"); // twords the sulfur
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[3]) == "c3"); // '                '
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[4]) == "ss"); // Linking sulfur
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[5]) == "c");  // Carboxylic acid
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[6]) == "o");  // Carboxylic acid
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[7]) == "o");  // Carboxylic acid
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[8]) == "c3"); // After sulfur
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[9]) == "c3"); // Ring carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[10])== "os"); // Ring oxygen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[11])== "c3"); // Ring carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[12])== "oh"); // hydroxyl
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[13])== "c3"); // Ring carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[14])== "oh"); // hydroxyl
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[15])== "c3"); // Final ring carbon
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[16])== "na"); // First aromatic nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[17])== "cc"); // Carbon in non-pure aro ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[18])== "nd"); // nitrogren in 5-membered ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[19])== "ca"); // Carbon in 6-membered aro ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[20])== "ca"); // Carbon in 6-membered aro ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[21])== "nh"); // Nitrogen next to ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[22])== "nb"); // 'Pure' aro nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[23])== "ca"); // Carbon in 6-membered aro ring
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[24])== "nb"); // 'Pure' aro nitrogen
    CHECK(atomtype_name_for_id<GAFF>(atomtypes[25])== "ca"); // Carbon in 6-membered aro ring
}
