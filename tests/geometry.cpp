#include "spear/Geometry.hpp"
#include "spear/Molecule.hpp"
#include <utility>
#include <numeric>

#include "chemfiles.hpp"

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

TEST_CASE("Kabsch") {
    Spear::Conformation in(100);
    Spear::Conformation out(100);
    Eigen::Quaternion<double> Q(1, 3, 5, 2);

    Q.normalize();
    auto R = Q.toRotationMatrix();
    double scale = 2.0;
    for (size_t row = 0; row < 3; row++) {
        for (size_t col = 0; col < in.size(); col++) {
            in[col][static_cast<int>(row)] = 
                std::sqrt(static_cast<double>(col)*1.0)/
                          (static_cast<double>(row) + 1.0) +
                std::log(2*static_cast<double>(row) + 10.0)/
                std::sqrt(1.0*static_cast<double>(col) + 4.0);
        }
    }

    Eigen::Vector3d S;
    S << -5, 6, -27;
    for (size_t col = 0; col < in.size(); col++) {
        out[col] = scale*R*in[col] + S;
    }

    Eigen::Affine3d A = Spear::kabsch(in, out);

    // See if we got the transform we expected
    CHECK((scale*R-A.linear()).cwiseAbs().maxCoeff() < 1e-13);
    CHECK((S-A.translation()).cwiseAbs().maxCoeff() < 1e-13);
}

TEST_CASE("RMSD") {
    chemfiles::Trajectory hiv1("data/1aaq_1.pdb");
    chemfiles::Trajectory hiv2("data/1aaq_2.pdb");
    Spear::Molecule mol2(hiv1.read());
    Spear::Molecule mol1(hiv2.read());

    double rmsd = Spear::rmsd(mol1.positions(), mol2.positions());
    CHECK(rmsd > 30); // should be about 30

    auto A = Spear::kabsch(mol1.positions(), mol2.positions());
    auto mol2_positions = mol2.positions();
    for (auto& p : mol2_positions) {
        p = A.linear() * (p - A.translation());
    }

    double rmsd2 = Spear::rmsd(mol1.positions(), mol2_positions);
    CHECK(rmsd2 < 1.0); // should be much smaller now
}
