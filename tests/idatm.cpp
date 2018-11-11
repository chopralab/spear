#include "spear/Molecule.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

auto idatm_name = atomtype_name_for_id<IDATM>;
auto idatm_type = atomtype_id_for_name<IDATM>;

TEST_CASE("IDATM") {
    SECTION("Atom Type Info") {
        CHECK(atomtype_id_count<IDATM>() == 150);
    }

    SECTION("Palmitic Acid") {
        auto traj = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol = Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

        auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
        CHECK(unique_types.size() == 3);

        CHECK(idatm_name(alltypes[16]) == "O2-");
        CHECK(idatm_name(alltypes[17]) == "O2-");
        CHECK(atomtype_name_for_id<IDATM>(alltypes[15]) == "Cac");

        // Check all but the last 3
        for (size_t i = 0; i < alltypes.size() - 3; ++i) {
            CHECK(atomtype_name_for_id<IDATM>(alltypes[i]) == "C3");
        }

        IDATM idatm2(mol);
        idatm2.type_atoms_order();
        auto alltypes_2 = idatm2.all_types();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Tibolone") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

        auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
        CHECK(unique_types.size() == 5);
        CHECK(unique_types.count(idatm_type("C3")) != 0);
        CHECK(unique_types.count(idatm_type("C2")) != 0);
        CHECK(unique_types.count(idatm_type("C1")) != 0);
        CHECK(unique_types.count(idatm_type("O2")) != 0);
        CHECK(unique_types.count(idatm_type("O3")) != 0);

        IDATM idatm2(mol);
        idatm2.type_atoms_order();
        auto alltypes_2 = idatm2.all_types();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Pazopanib") {
        auto traj = chemfiles::Trajectory("data/pazopanib.sdf");
        auto mol = Spear::Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

        auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
        CHECK(unique_types.size() == 6);
        CHECK(unique_types.count(idatm_type("C3")) != 0);
        CHECK(unique_types.count(idatm_type("Car")) != 0);
        CHECK(unique_types.count(idatm_type("N2")) != 0);
        CHECK(unique_types.count(idatm_type("Npl")) != 0);
        CHECK(unique_types.count(idatm_type("O3-")) != 0);
        CHECK(unique_types.count(idatm_type("Son")) != 0);

        auto rings = mol.rings();
        for (auto ring : rings) {
            for (auto atom : ring) {
                if (ring.size() == 6) {
                    if (mol[atom].atomic_number() == 6)
                        CHECK(idatm_name(alltypes[atom]) == "Car");
                    else if (mol[atom].atomic_number() == 7)
                        CHECK(idatm_name(alltypes[atom]) == "N2");
                    else
                        CHECK(0);
                } else if (ring.size() == 5) {
                    if (mol[atom].atomic_number() != 7) continue;
                    if (mol[atom].neighbor_count() == 3) {
                        CHECK(idatm_name(alltypes[atom]) == "Npl");
                    } else {
                        CHECK(idatm_name(alltypes[atom]) == "N2");
                    }
                }
            }
        }

        IDATM idatm2(mol);
        idatm2.type_atoms_order();
        auto alltypes_2 = idatm2.all_types();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Aromatics") {
        auto traj = chemfiles::Trajectory("data/aromatics.sdf");
        auto mol = Spear::Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

        auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
        CHECK(unique_types.size() == 10);
        CHECK(unique_types.count(idatm_type("Car")) != 0);
        CHECK(unique_types.count(idatm_type("Oar")) != 0);
        CHECK(unique_types.count(idatm_type("Oar+")) != 0);
        CHECK(unique_types.count(idatm_type("O3-")) != 0);
        CHECK(unique_types.count(idatm_type("N2")) != 0);
        CHECK(unique_types.count(idatm_type("N2+")) != 0);
        CHECK(unique_types.count(idatm_type("N1")) != 0);
        CHECK(unique_types.count(idatm_type("N1+")) != 0);
        CHECK(unique_types.count(idatm_type("Sar")) != 0);
        CHECK(unique_types.count(idatm_type("H")) != 0);

        IDATM idatm2(mol);
        idatm2.type_atoms_order();
        auto alltypes_2 = idatm2.all_types();
std::cout << "72: " << idatm_name(47) << " 73: " << idatm_name(48) << std::endl;
        CHECK(alltypes == alltypes_2);
	}

    SECTION("Oxides") {
        auto traj = chemfiles::Trajectory("data/oxides.sdf");
        auto mol = Spear::Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

		auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
		CHECK(unique_types.size() == 10);
        CHECK(unique_types.count(idatm_type("C3")) != 0);
        CHECK(unique_types.count(idatm_type("C2")) != 0);
        CHECK(unique_types.count(idatm_type("O3")) != 0);
        CHECK(unique_types.count(idatm_type("O3-")) != 0);
        CHECK(unique_types.count(idatm_type("Ng+")) != 0);
        CHECK(unique_types.count(idatm_type("Nox")) != 0);
        CHECK(unique_types.count(idatm_type("Pac")) != 0);
        CHECK(unique_types.count(idatm_type("Pox")) != 0);
        CHECK(unique_types.count(idatm_type("Sac")) != 0);
        CHECK(unique_types.count(idatm_type("Sxd")) != 0);

        IDATM idatm2(mol);
        idatm2.type_atoms_order();
        auto alltypes_2 = idatm2.all_types();
        CHECK(alltypes == alltypes_2);
	}

    SECTION("POB Problem case") {
        auto traj = chemfiles::Trajectory("data/pob.sdf");
        auto mol = Spear::Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

		auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
		CHECK(unique_types.size() == 8);
        CHECK(unique_types.count(idatm_type("O3")) != 0);
        CHECK(unique_types.count(idatm_type("O3-")) != 0);
        CHECK(unique_types.count(idatm_type("O2-")) != 0);
        CHECK(unique_types.count(idatm_type("Cac")) != 0);
        CHECK(unique_types.count(idatm_type("N3")) != 0);
        CHECK(unique_types.count(idatm_type("Pac")) != 0);
        CHECK(unique_types.count(idatm_type("Pox")) != 0);
        CHECK(unique_types.count(idatm_type("C3")) != 0);

        // The problem, no C2 should be presnt
        CHECK(unique_types.count(idatm_type("C2")) == 0);
    }

    SECTION("0T8 Problem case") {
        auto traj = chemfiles::Trajectory("data/0t8.sdf");
        auto mol = Spear::Molecule(traj.read());

        IDATM idatm(mol);
        idatm.type_atoms_3d();
        auto alltypes = idatm.all_types();
        CHECK(alltypes.size() == mol.size());

		auto unique_types = std::unordered_set<size_t>(alltypes.cbegin(), alltypes.cend());
		CHECK(unique_types.size() == 6);
        CHECK(unique_types.count(idatm_type("O2")) != 0);
        CHECK(unique_types.count(idatm_type("C2")) != 0);
        CHECK(unique_types.count(idatm_type("Npl")) != 0);
        CHECK(unique_types.count(idatm_type("Car")) != 0);
        CHECK(unique_types.count(idatm_type("N2")) != 0);
        CHECK(unique_types.count(idatm_type("C3")) != 0);

        auto rings = mol.rings();
        for (auto ring : rings) {
            for (auto atom : ring) {
                if (ring.size() == 5) {
                    CHECK(idatm_name(alltypes[atom]) != "Car");
                }
            }
        }
    }
}
