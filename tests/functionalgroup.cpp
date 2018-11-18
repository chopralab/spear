#include "spear/FunctionalGroup.hpp"
#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Functional Group") {
    SECTION("Tibolone") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::FunctionalGroup tib_OH("CCO");
        auto oh_grp = Spear::find_functional_groups(mol, tib_OH);
        CHECK(oh_grp.size() == 3);

        Spear::FunctionalGroup alkene("C=C");
        auto alkene_grp = Spear::find_functional_groups(mol, alkene);
        CHECK(alkene_grp.size() == 2); // two orders of the same bond

        Spear::FunctionalGroup ketone("CC(=O)C");
        auto ketone_grp = Spear::find_functional_groups(mol, ketone);
        CHECK(ketone_grp.size() == 2); // two orders of the same bond

        Spear::FunctionalGroup alkyne("CC#C"); // Force a order for the alkyne
        auto alkyne_grp = Spear::find_functional_groups(mol, alkyne);
        CHECK(alkyne_grp.size() == 1); // just one now

        // another example of multiple matches
        Spear::FunctionalGroup cyclo_pentane("C1CCCC1");
        auto cpentane_grp = Spear::find_functional_groups(mol, cyclo_pentane);
        CHECK(cpentane_grp.size() == 10); // Both ways... 5 times!

        Spear::FunctionalGroup cyclo_pentane_oh("OC1CCCC1");
        cpentane_grp = Spear::find_functional_groups(mol, cyclo_pentane_oh);
        CHECK(cpentane_grp.size() == 2); // Both ways, but not five ways!

        Spear::FunctionalGroup fused_ring("C12CCCC1CCCC2");
        auto fused_grp = Spear::find_functional_groups(mol, fused_ring);
        CHECK(fused_grp.size() == 2); // ring is topologically symmetrical

        Spear::FunctionalGroup fused_ring_oh("C12CCC(O)C1CCCC2");
        fused_grp = Spear::find_functional_groups(mol, fused_ring_oh);
        CHECK(fused_grp.size() == 1); // forced again
    }
}
