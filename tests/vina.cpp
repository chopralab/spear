#include "spear/Molecule.hpp"

#include "spear/scoringfunctions/VinaScore.hpp"
#include "spear/atomtypes/VinaType.hpp"
#include "spear/Grid.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

auto vina_name = atomtype_name_for_id<VinaType>;
auto vina_type = atomtype_id_for_name<VinaType>;

TEST_CASE("Protein-Ligand Score") {
    auto ptraj = chemfiles::Trajectory("data/3qox_pocket.pdb");
    auto protein = Molecule(ptraj.read());
    auto atomtype_name = protein.add_atomtype<VinaType>(VinaType::ACCEPTOR_WATER);

    auto ltraj = chemfiles::Trajectory("data/3qox_ligand.sdf");
    auto ligand = Molecule(ltraj.read());
    auto atomtype_name2 = ligand.add_atomtype<VinaType>();

    CHECK(atomtype_name == atomtype_name2);

    auto at = ligand.atomtype(atomtype_name);
    CHECK(at != nullptr);

    auto& atr = *at;
    CHECK(vina_name(atr[0]) == "N_D");
    CHECK(vina_name(atr[1]) == "C_P");
    CHECK(vina_name(atr[2]) == "C_H");
    CHECK(vina_name(atr[3]) == "C_P");
    CHECK(vina_name(atr[4]) == "S_P");
    CHECK(vina_name(atr[5]) == "C_P");
    CHECK(vina_name(atr[6]) == "O_A");
    CHECK(vina_name(atr[7]) == "O_A");
    CHECK(vina_name(atr[8]) == "C_P");
    CHECK(vina_name(atr[9]) == "C_P");
    CHECK(vina_name(atr[10])== "O_A");
    CHECK(vina_name(atr[11])== "C_P");
    CHECK(vina_name(atr[12])== "O_DA");
    CHECK(vina_name(atr[13])== "C_P");
    CHECK(vina_name(atr[14])== "O_DA");
    CHECK(vina_name(atr[15])== "C_P");
    CHECK(vina_name(atr[16])== "N_P");
    CHECK(vina_name(atr[17])== "C_P");
    CHECK(vina_name(atr[18])== "N_A");
    CHECK(vina_name(atr[19])== "C_P");
    CHECK(vina_name(atr[20])== "C_P");
    CHECK(vina_name(atr[21])== "N_D");
    CHECK(vina_name(atr[22])== "N_A");
    CHECK(vina_name(atr[23])== "C_P");
    CHECK(vina_name(atr[24])== "N_A");
    CHECK(vina_name(atr[25])== "C_P");
    CHECK(vina_name(atr[26])== "SKIP");

    auto grid = Grid(protein.positions());
    VinaScore scoring_func;
    auto thing = scoring_func.calculate_components(grid, protein, ligand);
    CHECK(std::fabs(thing.g1 - 0153.3658) < 1e-3);
    CHECK(std::fabs(thing.g2 - 1837.0294) < 1e-3);
    CHECK(std::fabs(thing.rep - 8.8736) < 1e-3);
    CHECK(std::fabs(thing.hydrogen - 12.9489) < 1e-3);
    CHECK(std::fabs(thing.hydrophobic - 1.7243) < 1e-3);
    CHECK(std::fabs(scoring_func.score(grid, protein, ligand) - -15.1394) < 1e-3);
}
