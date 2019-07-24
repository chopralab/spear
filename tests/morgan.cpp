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
    CHECK(fp_0.size() == 5);

    auto fp_1 = calcFingerprint(mol, 1);
    CHECK(fp_1.size() == 12);

    auto fp_2 = calcFingerprint(mol, 2);
    CHECK(fp_2.size() == 18);

    auto fp_3 = calcFingerprint(mol, 3);
    CHECK(fp_3.size() == 24);

    auto fp_4 = calcFingerprint(mol, 4);
    CHECK(fp_4.size() == 30);

    auto fp_5 = calcFingerprint(mol, 5);
    CHECK(fp_5.size() == 36);

    auto fp_6 = calcFingerprint(mol, 6);
    CHECK(fp_6.size() == 42);

    auto fp_7 = calcFingerprint(mol, 7);
    CHECK(fp_7.size() == 47);

    auto fp_8 = calcFingerprint(mol, 8);
    CHECK(fp_8.size() == 50);

    // Converged
    auto fp_9 = calcFingerprint(mol, 9);
    CHECK(fp_9.size() == 50);
}
