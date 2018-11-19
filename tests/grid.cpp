#include "spear/Grid.hpp"
#include "chemfiles.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

TEST_CASE("Grid Contrived") {
    auto points = std::vector<Vector3D> {
        {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}, 
        {2.5, 2.5, 2.5}, {3.0, 3.0, 3.0},
    };

    CHECK_THROWS(Grid(points, {1.0, 0.0, 0.0}));

    Grid grid(points);
    CHECK(grid.size() == 4);

    auto occupied = grid.occupied();
    CHECK(occupied.size() == 3);
    CHECK(occupied.count(GridPoint(0,0,0)));
    CHECK(occupied.count(GridPoint(1,1,1)));
    CHECK(occupied.count(GridPoint(2,2,2)));

    auto no_neighbors = grid.neighbors({1.0, 3.5, 0.0});
    CHECK(no_neighbors.size() == 0);

    auto has_neighbors1 = grid.neighbors({1.6, 2.1, 1.4}, 0.5);
    CHECK(has_neighbors1.size() == 1);

    auto has_neighbors3 = grid.neighbors({1.6, 2.1, 1.6}, 0.5);
    CHECK(has_neighbors3.size() == 3);
}

TEST_CASE("Grid protein") {
    auto prot = chemfiles::Trajectory("data/3qox_pocket.pdb");
    const auto protf = prot.read();

    Grid grid(protf.positions(), {2.0, 2.0, 2.0});
    CHECK(grid.size() == protf.size());
    auto occupied = grid.occupied();
    CHECK(occupied.size() < grid.size());
}
