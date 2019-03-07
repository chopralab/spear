#include "spear/FunctionalGroup.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Tibolone") {
    auto traj = chemfiles::Trajectory("data/tibolone.sdf");
    auto mol = Spear::Molecule(traj.read());

    Spear::FunctionalGroup tib_O("CCO");
    auto o_grp = Spear::find_functional_groups(mol, tib_O);
    CHECK(o_grp.size() == 5); // Bond order not set, matches C-O and C=O!!!

    Spear::FunctionalGroup tib_OH("C-C-O");
    auto oh_grp = Spear::find_functional_groups(mol, tib_OH);
    CHECK(oh_grp.size() == 3); // Only matches C-O!!!

    Spear::FunctionalGroup alkene("C=C");
    auto alkene_grp = Spear::find_functional_groups(mol, alkene);
    CHECK(alkene_grp.size() == 2); // two orders of the same bond

    Spear::FunctionalGroup ketone("C-C(=O)-C");
    auto ketone_grp = Spear::find_functional_groups(mol, ketone);
    CHECK(ketone_grp.size() == 2); // two orders of the same bond

    Spear::FunctionalGroup alkyne("C-C#C"); // Force a order for the alkyne
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

    Spear::FunctionalGroup wildcard("C-C(=*)-C");
    auto wild_grp = Spear::find_functional_groups(mol, wildcard);
    CHECK(wild_grp.size() == 6); // forced again
}

TEST_CASE("Pazopanib") {
    auto traj = chemfiles::Trajectory("data/pazopanib.sdf");
    auto mol = Spear::Molecule(traj.read());
    auto name = mol.add_atomtype<Spear::IDATM>(Spear::AtomType::GEOMETRY);
    mol.set_default_atomtype(name);

    SECTION("General Searching") {
        Spear::FunctionalGroup aromatic_C_N("C:N");
        auto a_C_N_grp = Spear::find_functional_groups(mol, aromatic_C_N);
        CHECK(a_C_N_grp.size() == 6); // 6 aromatic carbon-nitrogen bonds!

        Spear::FunctionalGroup aromatic_ring_sym("C1:N:C(-N):N:C:C:1");
        auto ring_sym_grp = Spear::find_functional_groups(mol, aromatic_ring_sym);
        CHECK(ring_sym_grp.size() == 2); // Symmetric, there's two possibilities

        Spear::FunctionalGroup aromatic_ring_sym_2("c1nc(-N)ncc1");
        auto ring_sym_grp_2 = Spear::find_functional_groups(mol, aromatic_ring_sym_2);
        CHECK(ring_sym_grp_2.size() == 2); // Statement should be the same as above.

        Spear::FunctionalGroup aromatic_ring_sym_broken("C1:N[A]:C(-N):N:C:C:1");
        auto ring_sym_grp_broken = Spear::find_functional_groups(mol, aromatic_ring_sym_broken);
        CHECK(ring_sym_grp_broken.size() == 0); // We forced something to be aliphatic

        Spear::FunctionalGroup aromatic_ring("C1:N:C(-N):N:C(-N):C:1");
        auto ring_grp = Spear::find_functional_groups(mol, aromatic_ring);
        CHECK(ring_grp.size() == 1); // No longer symmetric!

        Spear::FunctionalGroup sulfonamide("S(=O)(=O)-N");
        auto so2n_grp = Spear::find_functional_groups(mol, sulfonamide);
        CHECK(so2n_grp.size() == 2); // Oxygens are symmetric
    }

    SECTION("Explicit Bonds") {
        // Explicit Bonds
        Spear::FunctionalGroup N_with_3_bonds("[ND3]");
        auto N_3_grp = Spear::find_functional_groups(mol, N_with_3_bonds);
        CHECK(N_3_grp.size() == 2); // Only two nitrogens have three explicit(non-H) bonds

        Spear::FunctionalGroup N_with_3_bonds_aliphatic("[ND3A]");
        auto N_3_ali_grp = Spear::find_functional_groups(mol, N_with_3_bonds_aliphatic);
        CHECK(N_3_ali_grp.size() == 1); // But only one of those nitrogens is aliphatic!

        Spear::FunctionalGroup N_with_2_bonds("[ND2]");
        auto N_2_grp = Spear::find_functional_groups(mol, N_with_2_bonds);
        CHECK(N_2_grp.size() == 4); // Four nitrogens have two explicit(non-H) bonds

        Spear::FunctionalGroup N_with_2_bonds_aliphatic("[ND2A]");
        auto N_2_ali_grp = Spear::find_functional_groups(mol, N_with_2_bonds_aliphatic);
        CHECK(N_2_ali_grp.size() == 1); // But only one of those nitrogens is aliphatic!

        Spear::FunctionalGroup N_with_1_bonds("[ND1]");
        auto N_1_grp = Spear::find_functional_groups(mol, N_with_1_bonds);
        CHECK(N_1_grp.size() == 1); // sulfonamide nitrogen!
    }

    SECTION("Expected Bonds") {
        // Expected Bonds
        Spear::FunctionalGroup N_with_3_bonds("[NX3]");
        auto N_3_grp = Spear::find_functional_groups(mol, N_with_3_bonds);
        CHECK(N_3_grp.size() == 4); // Four nitrogens have 3 expected bonds

        Spear::FunctionalGroup N_with_2_bonds("[NX2]");
        auto N_2_grp = Spear::find_functional_groups(mol, N_with_2_bonds);
        CHECK(N_2_grp.size() == 3); // Three nitrogens have 2 expected bonds

        Spear::FunctionalGroup N_with_2_bonds_aliphatic("[NX2A]");
        auto N_2_ali_grp = Spear::find_functional_groups(mol, N_with_2_bonds_aliphatic);
        CHECK(N_2_ali_grp.size() == 0); // But none of them are aliphatic!

        Spear::FunctionalGroup N_with_2_bonds_aromatic("[NX2a]");
        auto N_2_aro_grp = Spear::find_functional_groups(mol, N_with_2_bonds_aromatic);
        CHECK(N_2_aro_grp.size() == 3); // They are aromatic!
    }
}
