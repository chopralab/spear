#include "spear/Octree.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using namespace Spear;

TEST_CASE("Octree Contrived") {
    auto points = Conformation {
        {-1.0,-1.0,-1.0}, {2.0, 2.0, 2.0}, 
        { 2.5, 2.5, 2.5}, {3.0, 3.0, 3.0},
    };

    auto bounds = Octree::bounds(points);
    CHECK(bounds.center[0] == 1.0);
    CHECK(bounds.center[1] == 1.0);
    CHECK(bounds.center[2] == 1.0);
    CHECK(bounds.radius == 4.0);

    Octree tree(points, 2, 10);
    CHECK(tree.size() == 4);

    CHECK(tree.neighbors({0.0, 0.0, 0.0}).size() == 1);
    CHECK(tree.neighbors({0.0, 0.0, 0.0})[0] == 0);

    CHECK(tree.neighbors({1.5, 1.5, 1.5}).size() == 1);
    CHECK(tree.neighbors({1.5, 1.5, 1.5})[0] == 1);

    CHECK(tree.neighbors({2.6, 2.6, 2.6}).size() == 2);
    CHECK(tree.neighbors({2.6, 2.6, 2.6})[0] == 2);
    CHECK(tree.neighbors({2.6, 2.6, 2.6})[1] == 3);
}
