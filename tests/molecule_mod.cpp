#include "spear/Molecule.hpp"
#include "spear/Geometry.hpp"

#ifndef M_PI
static const auto M_PI = std::acos(0.0) * 2;
#endif 


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
        auto traj = chemfiles::Trajectory("data/palmitic_2d.sdf");
        auto mol = Spear::Molecule(traj.read());

        CHECK(Spear::dimensionality(mol.positions()) == 2);

        auto n_atom = mol.add_atom_to(Spear::Element::N, 0);
        CHECK(mol[0].degree() == 2);
        auto mg_atom = mol.add_atom_to(Spear::Element::Mg, 0);
        CHECK(mol[0].degree() == 3);
        mol.add_atom_to(Spear::Element::H, 0);
        CHECK(mol[0].degree() == 4);

        auto& atom_types = *mol.get_default_atomtype();
        CHECK(atom_types.hybridization(n_atom) == Spear::Hybridization::SP3);
        CHECK(atom_types.hybridization(mg_atom)== Spear::Hybridization::UNKNOWN);
    }

    SECTION("Given Atom in 3D") {
        auto traj = chemfiles::Trajectory("data/tibolone.sdf");
        auto mol = Spear::Molecule(traj.read());

        CHECK(Spear::dimensionality(mol.positions()) == 3);

        // A methyl group
        auto h1 = mol.add_atom_to(Spear::Element::H, 0);
        auto h2 = mol.add_atom_to(Spear::Element::H, 0);
        auto h3 = mol.add_atom_to(Spear::Element::H, 0);
        CHECK(std::fabs(Spear::distance(mol[0].position(), h1.position()) - 1.06) < 1e-4);
        auto tetra_ang = std::fabs(Spear::angle(h2.position(), mol[0].position(), h1.position()));
        CHECK(std::fabs(tetra_ang * 180.0 / M_PI - 109.47) < 1e-1);
        auto dihedral_ang = std::fabs(Spear::dihedral(h2.position(), mol[0].position(), h1.position(), h3.position()));
        CHECK(std::fabs(dihedral_ang * 180.0 / M_PI - 120.0) < 1e-1);

        auto& atom_types = *mol.get_default_atomtype();
        CHECK(atom_types[h1] == Spear::Element::H);
        CHECK(atom_types[h2] == Spear::Element::H);
        CHECK(atom_types[h3] == Spear::Element::H);
        CHECK(atom_types[0]  == Spear::Element::C);

        CHECK(atom_types.hybridization(h1) == Spear::Hybridization::FORCED);
        CHECK(atom_types.hybridization(h2) == Spear::Hybridization::FORCED);
        CHECK(atom_types.hybridization(h3) == Spear::Hybridization::FORCED);

        auto h4 = mol.add_atom_to(Spear::Element::H, 21);
        CHECK(std::fabs(Spear::distance(mol[21].position(), h4.position()) - 0.99) < 1e-4);
        auto linear_ang = std::fabs(Spear::angle(h4.position(), mol[21].position(), mol[20].position()));
        CHECK(std::fabs(linear_ang * 180.0 / M_PI - 180.0) < 1e-1);

        // Full valence SP
        CHECK_THROWS(mol.add_atom_to(Spear::Element::H, 20));
        CHECK_THROWS(mol.add_atom_to(Spear::Element::H, 21));

        // Full valence SP2
        CHECK_THROWS(mol.add_atom_to(Spear::Element::H, 20));

        // Can't add a bond to an existing bond!
        CHECK_THROWS(mol.add_bond(20, 21, Spear::Bond::SINGLE));
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
        CHECK(mol.get_default_atomtype()->size() == mol.size());
    }

    SECTION("Add") {
        auto mol_3qox = Spear::Molecule(chemfiles::Trajectory("data/3qox_ligand.sdf").read());
        CHECK(mol_3qox.add_hydrogens() == 1); // TODO: Again a pH thing

        auto mol_tib = Spear::Molecule(chemfiles::Trajectory("data/tibolone.sdf").read());
        CHECK(mol_tib.add_hydrogens() == 28);

        auto mol_paz = Spear::Molecule(chemfiles::Trajectory("data/pazopanib.sdf").read());
        CHECK(mol_paz.add_hydrogens() == 23);
    }
}
