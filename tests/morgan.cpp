#include "spear/Molecule.hpp"
#include "spear/Fingerprint.hpp"
#include "spear/fingerprints/Morgan.hpp"

#include "chemfiles.hpp"
#include <utility>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Palmitic acid") {
    auto traj = chemfiles::Trajectory("data/palmitic.sdf");
    auto mol = Spear::Molecule(traj.read());

    auto fp_0 = calcFingerprint(mol, 0);
    CHECK(fp_0.non_zero_values() == 5);

    auto fp_1 = calcFingerprint(mol, 1);
    CHECK(fp_1.non_zero_values() == 12);

    auto fp_2 = calcFingerprint(mol, 2);
    CHECK(fp_2.non_zero_values() == 18);

    auto fp_3 = calcFingerprint(mol, 3);
    CHECK(fp_3.non_zero_values() == 24);

    auto fp_4 = calcFingerprint(mol, 4);
    CHECK(fp_4.non_zero_values() == 30);

    auto fp_5 = calcFingerprint(mol, 5);
    CHECK(fp_5.non_zero_values() == 36);

    auto fp_6 = calcFingerprint(mol, 6);
    CHECK(fp_6.non_zero_values() == 42);

    auto fp_7 = calcFingerprint(mol, 7);
    CHECK(fp_7.non_zero_values() == 47);

    auto fp_8 = calcFingerprint(mol, 8);
    CHECK(fp_8.non_zero_values() == 50);

    // Converged
    auto fp_9 = calcFingerprint(mol, 9);
    CHECK(fp_9.non_zero_values() == 50);
}
