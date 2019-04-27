#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"

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
	CHECK(sybyl.size() == mol.size());

	auto unique_types = std::unordered_set<size_t>(sybyl.cbegin(), sybyl.cend());
    CHECK(unique_types.size() == 12);

    CHECK(sybyl_name(sybyl[0]) == "C.ar");
    CHECK(sybyl.is_aromatic(0));
    CHECK(sybyl.hybridization(0) == Hybridization::SP2);
    CHECK(sybyl_name(sybyl[1]) == "C.ar");
    CHECK(sybyl_name(sybyl[2]) == "N.ar");
    CHECK(sybyl.is_aromatic(2));
    CHECK(sybyl.hybridization(2) == Hybridization::SP2);
    CHECK(sybyl_name(sybyl[3]) == "C.ar");
    CHECK(sybyl_name(sybyl[4]) == "C.ar");
    CHECK(sybyl_name(sybyl[5]) == "C.ar");
    CHECK(sybyl_name(sybyl[6]) == "N.pl3");
    CHECK(!sybyl.is_aromatic(6));
    CHECK(sybyl.hybridization(6) == Hybridization::SP3);
    CHECK(sybyl_name(sybyl[7]) == "C.3");
    CHECK(sybyl_name(sybyl[8]) == "N.pl3");
    CHECK(sybyl_name(sybyl[9]) == "C.2");
    CHECK(sybyl_name(sybyl[10]) == "N.am");
    CHECK(sybyl_name(sybyl[11]) == "O.2");
    CHECK(sybyl_name(sybyl[12]) == "C.cat");
    CHECK(sybyl.hybridization(12) == Hybridization::SP2);
    CHECK(sybyl_name(sybyl[13]) == "N.pl3");
    CHECK(sybyl_name(sybyl[14]) == "N.pl3");
    CHECK(sybyl_name(sybyl[15]) == "C.2");
    CHECK(sybyl_name(sybyl[16]) == "C.2");
    CHECK(sybyl_name(sybyl[17]) == "C.3");
    CHECK(sybyl_name(sybyl[18]) == "C.3");
    CHECK(sybyl_name(sybyl[19]) == "C.1");
    CHECK(sybyl.hybridization(19) == Hybridization::SP);
    CHECK(sybyl_name(sybyl[20]) == "N.1");
    CHECK(sybyl_name(sybyl[21]) == "C.3");
    CHECK(sybyl_name(sybyl[22]) == "Fe");
    CHECK(sybyl.hybridization(22) == Hybridization::FORCED);
    CHECK(sybyl_name(sybyl[23]) == "C.3");
    CHECK(sybyl_name(sybyl[24]) == "H");
}
