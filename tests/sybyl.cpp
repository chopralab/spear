#include "spear/Molecule.hpp"
#include "spear/atomtypes/Sybyl.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

auto sybyl_name = atomtype_name_for_id<Sybyl>;
auto sybyl_type = atomtype_id_for_name<Sybyl>;

TEST_CASE("Sybyl count") {
    CHECK(atomtype_id_count<Sybyl>() == 128);
}

TEST_CASE("SYBYL") {
    auto traj = chemfiles::Trajectory("data/sybyl_test.mol2");
    auto mol = Molecule(traj.read());

    Sybyl sybyl(mol, AtomType::TOPOLOGY);
	auto alltypes = sybyl.all_types();
	CHECK(alltypes.size() == mol.size());

	auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
    CHECK(unique_types.size() == 12);

    CHECK(sybyl_name(alltypes[0]) == "C.ar");
    CHECK(sybyl_name(alltypes[1]) == "C.ar");
    CHECK(sybyl_name(alltypes[2]) == "N.ar");
    CHECK(sybyl_name(alltypes[3]) == "C.ar");
    CHECK(sybyl_name(alltypes[4]) == "C.ar");
    CHECK(sybyl_name(alltypes[5]) == "C.ar");
    CHECK(sybyl_name(alltypes[6]) == "N.pl3");
    CHECK(sybyl_name(alltypes[7]) == "C.3");
    CHECK(sybyl_name(alltypes[8]) == "N.pl3");
    CHECK(sybyl_name(alltypes[9]) == "C.2");
    CHECK(sybyl_name(alltypes[10]) == "N.am");
    CHECK(sybyl_name(alltypes[11]) == "O.2");
    CHECK(sybyl_name(alltypes[12]) == "C.cat");
    CHECK(sybyl_name(alltypes[13]) == "N.pl3");
    CHECK(sybyl_name(alltypes[14]) == "N.pl3");
    CHECK(sybyl_name(alltypes[15]) == "C.2");
    CHECK(sybyl_name(alltypes[16]) == "C.2");
    CHECK(sybyl_name(alltypes[17]) == "C.3");
    CHECK(sybyl_name(alltypes[18]) == "C.3");
    CHECK(sybyl_name(alltypes[19]) == "C.1");
    CHECK(sybyl_name(alltypes[20]) == "N.1");
    CHECK(sybyl_name(alltypes[21]) == "C.3");
    CHECK(sybyl_name(alltypes[22]) == "Fe");
    CHECK(sybyl_name(alltypes[23]) == "C.3");
    CHECK(sybyl_name(alltypes[24]) == "H");
}
