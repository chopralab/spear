#include "spear/partialcharges/Gasteiger.hpp"
#include "spear/Molecule.hpp"

#include "chemfiles.hpp"

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

TEST_CASE("Methane") {
    Spear::Molecule mol;
    mol.add_atom(Spear::Element::C, {0.00, 0.00, 0.00});
    mol.add_hydrogens();
    CHECK(mol.size() == 5);
    CHECK(mol[0].hybridization() == Spear::Hybridization::SP3);

    Spear::GasteigerCharge gc(mol);
    CHECK(gc.name() == "gasteiger");
    CHECK(std::abs(gc[0] - -0.0780) < 1e-3);
    CHECK(std::abs(gc[1] -  0.0195) < 1e-3);
    CHECK(std::abs(gc[2] -  0.0195) < 1e-3);
    CHECK(std::abs(gc[3] -  0.0195) < 1e-3);
    CHECK(std::abs(gc[4] -  0.0195) < 1e-3);
}

TEST_CASE("Fluoromethane") {
    Spear::Molecule mol;
    mol.add_atom(Spear::Element::C, {0.00, 0.00, 0.00});
    mol.add_atom_to(Spear::Element::F, 0);
    mol.add_hydrogens();
    CHECK(mol.size() == 5);

    Spear::GasteigerCharge gc(mol);
    CHECK(std::abs(gc[0] -  0.079) < 1e-3);
    CHECK(std::abs(gc[1] - -0.253) < 1e-3);
    CHECK(std::abs(gc[2] -  0.058) < 1e-3);
    CHECK(std::abs(gc[3] -  0.058) < 1e-3);
    CHECK(std::abs(gc[4] -  0.058) < 1e-3);
}

TEST_CASE("Ethane") {
    Spear::Molecule mol;
    mol.add_atom(Spear::Element::C, {0.00, 0.00, 0.00});
    mol.add_atom_to(Spear::Element::C, 0);
    mol.add_hydrogens();
    CHECK(mol.size() == 8);

    Spear::GasteigerCharge gc(mol);
    CHECK(std::abs(gc[0] - -0.0680) < 1e-3);
    CHECK(std::abs(gc[1] - -0.0680) < 1e-3);
    CHECK(std::abs(gc[2] -  0.0226) < 1e-3);
    CHECK(std::abs(gc[3] -  0.0226) < 1e-3);
    CHECK(std::abs(gc[4] -  0.0226) < 1e-3);
    CHECK(std::abs(gc[5] -  0.0226) < 1e-3);
    CHECK(std::abs(gc[6] -  0.0226) < 1e-3);
    CHECK(std::abs(gc[7] -  0.0226) < 1e-3);
}

TEST_CASE("Fluoroethan") {
    Spear::Molecule mol;
    mol.add_atom(Spear::Element::C, {0.00, 0.00, 0.00});
    mol.add_atom_to(Spear::Element::C, 0);
    mol.add_atom_to(Spear::Element::F, 0);
    mol.add_hydrogens();
    CHECK(mol.size() == 8);

    Spear::GasteigerCharge gc(mol);
    CHECK(std::abs(gc[0] -  0.0870) < 1e-3);
    CHECK(std::abs(gc[1] - -0.0370) < 1e-3);
    CHECK(std::abs(gc[2] - -0.2493) < 1e-3);
    CHECK(std::abs(gc[3] -  0.0611) < 1e-3);
    CHECK(std::abs(gc[4] -  0.0611) < 1e-3);
    CHECK(std::abs(gc[5] -  0.0254) < 1e-3);
    CHECK(std::abs(gc[6] -  0.0254) < 1e-3);
    CHECK(std::abs(gc[7] -  0.0254) < 1e-3);
}

TEST_CASE("Tibolone") {
    auto traj = chemfiles::Trajectory("data/tibolone.sdf");
    auto mol = Spear::Molecule(traj.read());
    Spear::GasteigerCharge gc(mol);

    CHECK(std::abs(gc[0] -  0.0003) < 1e-3);
    CHECK(std::abs(gc[1] -  0.0043) < 1e-3);
    CHECK(std::abs(gc[2] -  0.0218) < 1e-3);
    CHECK(std::abs(gc[3] - -0.0448) < 1e-3);
    CHECK(std::abs(gc[4] - -0.0516) < 1e-3);
    CHECK(std::abs(gc[5] -  0.0291) < 1e-3);
    CHECK(std::abs(gc[6] -  0.0692) < 1e-3);
    CHECK(std::abs(gc[7] -  0.1550) < 1e-3);
    CHECK(std::abs(gc[8] - -0.2957) < 1e-3);
    CHECK(std::abs(gc[9] -  0.0859) < 1e-3);
    CHECK(std::abs(gc[10]-  0.0177) < 1e-3);
    CHECK(std::abs(gc[11]-  0.0041) < 1e-3);
    CHECK(std::abs(gc[12]-  0.0037) < 1e-3);
    CHECK(std::abs(gc[13]-  0.0366) < 1e-3);
    CHECK(std::abs(gc[14]-  0.0034) < 1e-3);
    CHECK(std::abs(gc[15]-  0.0036) < 1e-3);
    CHECK(std::abs(gc[16]-  0.0037) < 1e-3);
    CHECK(std::abs(gc[17]-  0.0430) < 1e-3);
    CHECK(std::abs(gc[18]-  0.1797) < 1e-3);
    CHECK(std::abs(gc[19]- -0.2124) < 1e-3);
    CHECK(std::abs(gc[20]- -0.0491) < 1e-3);
    CHECK(std::abs(gc[21]- -0.0117) < 1e-3);
    CHECK(std::abs(gc[22]-  0.0041) < 1e-3);
}
