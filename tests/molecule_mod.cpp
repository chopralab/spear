#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Graph_impl.hpp"
#include "spear/Geometry.hpp"

#include "chemfiles.hpp"
#include <utility>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Add Atoms and Bonds") {
    SECTION("Given Point") {
        auto traj = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol = Spear::Molecule(traj.read());

        auto desc = mol.add_atom(Spear::Element::Cl, {0.7011, -5.3628, 8.4904});
        CHECK(mol.size() == 19);
        CHECK(desc.atomic_number() == Spear::Element::Cl);
        CHECK(desc == mol[18]);
        CHECK(desc.position() == Eigen::Vector3d{0.7011, -5.3628, 8.4904});
    }

    SECTION("Given Atom in 2D") {
        auto traj = chemfiles::Trajectory("data/palmitic.sdf");
        auto mol = Spear::Molecule(traj.read());

        CHECK(mol.dimensionality() == 2);

        mol.add_atom_to(Spear::Element::N, 0);
        CHECK(mol[0].neighbor_count() == 2);
        mol.add_atom_to(Spear::Element::H, 0);
        CHECK(mol[0].neighbor_count() == 3);
        mol.add_atom_to(Spear::Element::H, 0);
        CHECK(mol[0].neighbor_count() == 4);
    }

    SECTION("Given Atom in 3D") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        CHECK(mol.dimensionality() == 3);

        // A methyl group
        auto h1 = mol.add_atom_to(Spear::Element::H, 0);
        auto h2 = mol.add_atom_to(Spear::Element::H, 0);
        auto h3 = mol.add_atom_to(Spear::Element::H, 0);
        CHECK(std::fabs(Spear::distance(mol[0].position(), h1.position()) - 0.99) < 1e-4);
        auto tetra_ang = std::fabs(Spear::angle(h2.position(), mol[0].position(), h1.position()));
        CHECK(std::fabs(tetra_ang * 180.0 / M_PI - 109.47) < 1e-1);
        auto dihedral_ang = std::fabs(Spear::dihedral(h2.position(), mol[0].position(), h1.position(), h3.position()));
        CHECK(std::fabs(dihedral_ang * 180.0 / M_PI - 120.0) < 1e-1);

        auto h4 = mol.add_atom_to(Spear::Element::H, 21);
        CHECK(std::fabs(Spear::distance(mol[21].position(), h4.position()) - 0.99) < 1e-4);
        auto linear_ang = std::fabs(Spear::angle(h4.position(), mol[21].position(), mol[20].position()));
        CHECK(std::fabs(linear_ang * 180.0 / M_PI - 180.0) < 1e-1);
    }
}

TEST_CASE("Hydrogens") {
    SECTION("Remove") {
        auto traj = chemfiles::Trajectory("data/3qox_ligand.sdf");
        auto mol = Spear::Molecule(traj.read());

        size_t explicit_hs = 0, actual_hs = 0, implicit_hs = 0;
        for (auto av : mol) {
            explicit_hs += av.explicit_hydrogens();
            implicit_hs += av.implicit_hydrogens();
            if (av.atomic_number() == 1) ++actual_hs;
        }
        CHECK(actual_hs == 20);
        CHECK(implicit_hs == 1); // TODO: Fix carbonyl!
        CHECK(explicit_hs == actual_hs);

        mol.remove_hydrogens();

        explicit_hs = 0, actual_hs = 0, implicit_hs = 0;
        for (auto av : mol) {
            explicit_hs += av.explicit_hydrogens();
            implicit_hs += av.implicit_hydrogens();
            if (av.atomic_number() == 1) ++actual_hs;
        }
        CHECK(actual_hs == 0);
        CHECK(explicit_hs == actual_hs);
        CHECK(implicit_hs == 20);
        CHECK(mol.size() == 26);
    }
}
