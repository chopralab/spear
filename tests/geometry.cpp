#include "spear/Geometry.hpp"
#include <utility>
#include <numeric>

#define CATCH_CONFIG_MAIN
#include <catch.hpp>

using Spear::dimensionality;

TEST_CASE("Dimensionality") {
    SECTION("Point") {
        std::vector<Spear::Vector3d> point{ {1, 2, 3} };
        CHECK(dimensionality(point) == 0);

        point.emplace_back(1, 2, 3);
        CHECK(dimensionality(point) == 0);

        point.emplace_back(1, 2, 3);
        CHECK(dimensionality(point) == 0);

        point.emplace_back(1, 2, 3);
        CHECK(dimensionality(point) == 0);

        point.emplace_back(2, 3, 4);
        CHECK(dimensionality(point) == 1); // No longer 0 dimensions

        point.emplace_back(3, 4, 5);
        CHECK(dimensionality(point) == 1);

        point.emplace_back(0, 0, 0);
        CHECK(dimensionality(point) == 2); // No longer linear

        point.emplace_back(8, 8, 7);
        CHECK(dimensionality(point) == 3); // No longer planar
    }

    SECTION("Linear") {
        std::vector<Spear::Vector3d> line{ {1, 2, 3},
                                           {2, 3, 4}};

        CHECK(dimensionality(line) == 1);

        line.emplace_back(3, 4, 5);
        CHECK(dimensionality(line) == 1);

        line.emplace_back(4, 5, 6);
        CHECK(dimensionality(line) == 1);

        line.emplace_back(6, 8, 8);
        CHECK(dimensionality(line) == 2); // No longer linear

        line.emplace_back(6, 7, 8);
        CHECK(dimensionality(line) == 2); // Still planar

        line.emplace_back(7, 8, 8);
        CHECK(dimensionality(line) == 3); // No longer planar
    }

    SECTION("Planar") {
        std::vector<Spear::Vector3d> plane{ {0, 0, 0},
                                            {1, 0, 0},
                                            {0, 1, 0}};

        // L-shaped
        CHECK(dimensionality(plane) == 2);

        // T-shaped
        plane.emplace_back(-1, 0, 0);
        CHECK(dimensionality(plane) == 2);

        // No longer planer
        plane.emplace_back(0, 0, 1);
        CHECK(dimensionality(plane) == 3);        
    }
}
