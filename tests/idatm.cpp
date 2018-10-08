#include "spear/Molecule.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("IDATM") {
    SECTION("Palmitic Acid") {
        auto traj = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 3);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());

        CHECK(idatm.name(alltypes[16]) == "O2-");
        CHECK(idatm.name(alltypes[17]) == "O2-");
        CHECK(idatm.name(alltypes[15]) == "Cac");

        // Check all but the last 3
        for (size_t i = 0; i < alltypes.size() - 3; ++i) {
            CHECK(idatm.name(alltypes[i]) == "C3");
        }

        Spear::IDATM idatm_bo;
        idatm_bo.type_atoms_order(mol);
        auto alltypes_2 = idatm_bo.all_ids();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Tibolone") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 5);
        CHECK(unique_types.count(idatm.id("C3")) != 0);
        CHECK(unique_types.count(idatm.id("C2")) != 0);
        CHECK(unique_types.count(idatm.id("C1")) != 0);
        CHECK(unique_types.count(idatm.id("O2")) != 0);
        CHECK(unique_types.count(idatm.id("O3")) != 0);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());

        Spear::IDATM idatm_bo;
        idatm_bo.type_atoms_order(mol);
        auto alltypes_2 = idatm_bo.all_ids();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Pazopanib") {
        auto traj = chemfiles::Trajectory("data/pazopanib.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 6);
        CHECK(unique_types.count(idatm.id("C3")) != 0);
        CHECK(unique_types.count(idatm.id("Car")) != 0);
        CHECK(unique_types.count(idatm.id("N2")) != 0);
        CHECK(unique_types.count(idatm.id("Npl")) != 0);
        CHECK(unique_types.count(idatm.id("O3-")) != 0);
        CHECK(unique_types.count(idatm.id("Son")) != 0);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());

        auto rings = mol.rings();
        for (auto ring : rings) {
            for (auto atom : ring) {
                if (ring.size() == 6) {
                    if (mol[atom].atomic_number() == 6)
                        CHECK(idatm.name(alltypes[atom]) == "Car");
                    else if (mol[atom].atomic_number() == 7)
                        CHECK(idatm.name(alltypes[atom]) == "N2");
                    else
                        CHECK(0);
                } else if (ring.size() == 5) {
                    if (mol[atom].atomic_number() != 7) continue;
                    if (mol[atom].neighbor_count() == 3) {
                        CHECK(idatm.name(alltypes[atom]) == "Npl");
                    } else {
                        CHECK(idatm.name(alltypes[atom]) == "N2");
                    }
                }
            }
        }

        Spear::IDATM idatm_bo;
        idatm_bo.type_atoms_order(mol);
        auto alltypes_2 = idatm_bo.all_ids();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Aromatics") {
        auto traj = chemfiles::Trajectory("data/aromatics.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 10);
        CHECK(unique_types.count(idatm.id("Car")) != 0);
        CHECK(unique_types.count(idatm.id("Oar")) != 0);
        CHECK(unique_types.count(idatm.id("Oar+")) != 0);
        CHECK(unique_types.count(idatm.id("O3-")) != 0);
        CHECK(unique_types.count(idatm.id("N2")) != 0);
        CHECK(unique_types.count(idatm.id("N2+")) != 0);
        CHECK(unique_types.count(idatm.id("N1")) != 0);
        CHECK(unique_types.count(idatm.id("N1+")) != 0);
        CHECK(unique_types.count(idatm.id("Sar")) != 0);
        CHECK(unique_types.count(idatm.id("H")) != 0);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());

        Spear::IDATM idatm_bo;
        idatm_bo.type_atoms_order(mol);
        auto alltypes_2 = idatm_bo.all_ids();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("Oxides") {
        auto traj = chemfiles::Trajectory("data/oxides.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 10);
        CHECK(unique_types.count(idatm.id("C3")) != 0);
        CHECK(unique_types.count(idatm.id("C2")) != 0);
        CHECK(unique_types.count(idatm.id("O3")) != 0);
        CHECK(unique_types.count(idatm.id("O3-")) != 0);
        CHECK(unique_types.count(idatm.id("Ng+")) != 0);
        CHECK(unique_types.count(idatm.id("Nox")) != 0);
        CHECK(unique_types.count(idatm.id("Pac")) != 0);
        CHECK(unique_types.count(idatm.id("Pox")) != 0);
        CHECK(unique_types.count(idatm.id("Sac")) != 0);
        CHECK(unique_types.count(idatm.id("Sxd")) != 0);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());

        Spear::IDATM idatm_bo;
        idatm_bo.type_atoms_order(mol);
        auto alltypes_2 = idatm_bo.all_ids();
        CHECK(alltypes == alltypes_2);
    }

    SECTION("POB Problem case") {
        auto traj = chemfiles::Trajectory("data/pob.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 8);
        CHECK(unique_types.count(idatm.id("O3")) != 0);
        CHECK(unique_types.count(idatm.id("O3-")) != 0);
        CHECK(unique_types.count(idatm.id("O2-")) != 0);
        CHECK(unique_types.count(idatm.id("Cac")) != 0);
        CHECK(unique_types.count(idatm.id("N3")) != 0);
        CHECK(unique_types.count(idatm.id("Pac")) != 0);
        CHECK(unique_types.count(idatm.id("Pox")) != 0);
        CHECK(unique_types.count(idatm.id("C3")) != 0);

        // The problem, no C2 should be presnt
        CHECK(unique_types.count(idatm.id("C2")) == 0);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());
    }

    SECTION("0T8 Problem case") {
        auto traj = chemfiles::Trajectory("data/0t8.sdf");
        auto mol = Spear::Molecule(traj.read());

        Spear::IDATM idatm;
        idatm.type_atoms_3d(mol);

        auto unique_types = idatm.unique_ids();
        CHECK(unique_types.size() == 6);
        CHECK(unique_types.count(idatm.id("O2")) != 0);
        CHECK(unique_types.count(idatm.id("C2")) != 0);
        CHECK(unique_types.count(idatm.id("Npl")) != 0);
        CHECK(unique_types.count(idatm.id("Car")) != 0);
        CHECK(unique_types.count(idatm.id("N2")) != 0);
        CHECK(unique_types.count(idatm.id("C3")) != 0);

        auto alltypes = idatm.all_ids();
        CHECK(alltypes.size() == mol.size());

        auto rings = mol.rings();
        for (auto ring : rings) {
            for (auto atom : ring) {
                if (ring.size() == 5) {
                    CHECK(idatm.name(alltypes[atom]) != "Car");
                }
            }
        }
    }
}
