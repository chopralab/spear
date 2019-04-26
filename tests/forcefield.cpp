#include "spear/Forcefield.hpp"

#include <unordered_set>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Bonds") {
    std::unordered_set<std::array<size_t, 2>,
                       Spear::Forcefield::bond_type_hash,
                       Spear::Forcefield::bond_type_equal> bond_set;

    bond_set.insert({0, 1});
    bond_set.insert({1, 0});
    CHECK(bond_set.size() == 1);
}

TEST_CASE("Angles") {
    std::unordered_set<std::array<size_t, 3>,
                       Spear::Forcefield::angle_type_hash,
                       Spear::Forcefield::angle_type_equal> angle_set;

    angle_set.insert({0, 1, 2});
    angle_set.insert({2, 1, 0});
    CHECK(angle_set.size() == 1);
}

TEST_CASE("Torsions") {
    std::unordered_set<std::array<size_t, 4>,
                       Spear::Forcefield::torsion_type_hash,
                       Spear::Forcefield::torsion_type_equal> proper_set;

    Spear::Forcefield::torsion_type_hash hasher;
    Spear::Forcefield::torsion_type_equal checker;

    CHECK(hasher({3, 1, 2, 0}) == hasher({3, 1, 2, 0}));
    CHECK(hasher({3, 1, 2, 0}) == hasher({0, 2, 1, 3}));
    CHECK(hasher({0, 1, 2, 0}) == hasher({0, 2, 1, 0}));
    CHECK(hasher({0, 1, 1, 0}) == hasher({0, 1, 1, 0}));

    CHECK(hasher({0, 1, 1, 0}) != hasher({3, 1, 1, 3}));
    CHECK(hasher({0, 1, 1, 0}) != hasher({0, 2, 2, 0}));

    CHECK(checker({3, 1, 2, 0}, {0, 2, 1, 3}));
    CHECK(checker({0, 1, 2, 0}, {0, 2, 1, 0}));
    CHECK(checker({0, 2, 1, 0}, {0, 1, 2, 0}));

    proper_set.insert({0, 1, 2, 3});
    proper_set.insert({3, 2, 1, 0});
    CHECK(proper_set.size() == 1);

    proper_set.insert({3, 1, 2, 0});
    proper_set.insert({0, 2, 1, 3});
    CHECK(proper_set.size() == 2);

    proper_set.insert({0, 1, 1, 0});
    proper_set.insert({1, 0, 0, 1});
    CHECK(proper_set.size() == 4);
}
