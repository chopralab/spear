#include "spear/Molecule.hpp"
#include "spear/atomtypes/Sybyl.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("SYBYL") {
    auto traj = chemfiles::Trajectory("data/sybyl_test.mol2");
    auto mol = Spear::Molecule(traj.read());

    Spear::Sybyl sybyl;
    sybyl.type_atoms_order(mol);

    auto unique_types = sybyl.unique_ids();
    CHECK(unique_types.size() == 12);

    auto alltypes = sybyl.all_ids();
    CHECK(alltypes.size() == mol.size());

    CHECK(sybyl.name(alltypes[0]) == "C.ar");
    CHECK(sybyl.name(alltypes[1]) == "C.ar");
    CHECK(sybyl.name(alltypes[2]) == "N.ar");
    CHECK(sybyl.name(alltypes[3]) == "C.ar");
    CHECK(sybyl.name(alltypes[4]) == "C.ar");
    CHECK(sybyl.name(alltypes[5]) == "C.ar");
    CHECK(sybyl.name(alltypes[6]) == "N.pl3");
    CHECK(sybyl.name(alltypes[7]) == "C.3");
    CHECK(sybyl.name(alltypes[8]) == "N.pl3");
    CHECK(sybyl.name(alltypes[9]) == "C.2");
    CHECK(sybyl.name(alltypes[10]) == "N.am");
    CHECK(sybyl.name(alltypes[11]) == "O.2");
    CHECK(sybyl.name(alltypes[12]) == "C.cat");
    CHECK(sybyl.name(alltypes[13]) == "N.pl3");
    CHECK(sybyl.name(alltypes[14]) == "N.pl3");
    CHECK(sybyl.name(alltypes[15]) == "C.2");
    CHECK(sybyl.name(alltypes[16]) == "C.2");
    CHECK(sybyl.name(alltypes[17]) == "C.3");
    CHECK(sybyl.name(alltypes[18]) == "C.3");
    CHECK(sybyl.name(alltypes[19]) == "C.1");
    CHECK(sybyl.name(alltypes[20]) == "N.1");
    CHECK(sybyl.name(alltypes[21]) == "C.3");
    CHECK(sybyl.name(alltypes[22]) == "Fe");
    CHECK(sybyl.name(alltypes[23]) == "C.3");
    CHECK(sybyl.name(alltypes[24]) == "H");
}
