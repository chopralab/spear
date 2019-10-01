// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include "spear/Geometry.hpp"

using namespace Spear;

static double cos_diff(const Vector3d& u, const Vector3d& v) {
    return 1.0 - std::abs(u.dot(v) / (u.norm() * v.norm()));
}

size_t Spear::dimensionality(const std::vector<Vector3d>& positions, double eps) {
    
    auto size = positions.size();
    if (size == 0 || size == 1) {
        return 0;
    }

    // All the points could still be the same...
    bool is_zero_dimensions = distance(positions[0], positions[1]) < eps;

    // Any two points are on the same line
    if (size == 2) {
        return static_cast<size_t>(!is_zero_dimensions);
    }

    if (size == 3) {
        if (is_zero_dimensions && distance(positions[0], positions[1]) < eps) {
            return 0;
        }

        Vector3d lin_vec = positions[1] - positions[0];
        Vector3d plane_vec = positions[2] - positions[0];
        if (cos_diff(lin_vec, plane_vec) < eps) {
            return 1; // It's still linear!
        } else {
            return 2; // Any three points makes a plane
        }
    }

    // These two points, along with the first point, must be nonlinear
    // this allows one to calculate the nonplanar values one it is shown the
    // configuration is non-linear
    size_t nonlin1 = 0, nonlin2 = 0;

    bool is_linear = true;
    for (size_t i = 3; i < size; ++i) {
        if (is_zero_dimensions) {
            is_zero_dimensions = distance(positions[0], positions[i]) < eps;
            continue;
        }

        Vector3d curr_vec = positions[i] - positions[0];
        if (is_linear) {
            // Check all previous points for non-linearity
            // If we find something non-linear, we store this point
            // and the other non-linear point to check for non-planarity
            for (size_t j = 0; j < i; ++j) {
                if (cos_diff(positions[j] - positions[0], curr_vec) > eps) {
                    nonlin1 = j;
                    nonlin2 = i;
                    is_linear = false;
                    break;
                }
            }

            continue;
        }

        auto nplane = nonplanar(positions[0], positions[i],
                                positions[nonlin1], positions[nonlin2]
        );
        if (std::abs(nplane) > eps) {
            return 3;
        }
    }

    if (is_zero_dimensions) {
        return 0;
    } else if (is_linear) {
        return 1;
    } else {
        return 2;
    }
}


Eigen::Affine3d Spear::kabsch(const Conformation& conform1,
                              const Conformation& conform2) {

    // Default output
    Eigen::Affine3d A;
    A.linear() = Eigen::Matrix3d::Identity(3, 3);
    A.translation() = Eigen::Vector3d::Zero();

    if (conform1.size() != conform2.size()) {
        throw std::runtime_error("kabsch(): input data mis-match");
    }

    // First find the scale, by finding the ratio of sums of some distances,
    // then bring the datasets to the same scale.
    double dist1 = 0, dist2 = 0;
    for (size_t i = 0; i < conform1.size()-1; i++) {
        dist1 += (conform1[i+1] - conform1[i]).norm();
        dist2 += (conform2[i+1] - conform2[i]).norm();
    }
    if (dist1 <= 0 || dist2 <= 0) {
        return A;
    }
    double scale = dist2/dist1;

    // Find the centroids then shift to the origin
    Eigen::Vector3d ctr1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d ctr2 = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < conform1.size(); i++) {
        ctr1 += conform1[i];
        ctr2 += conform2[i];
    }
    ctr1 /= static_cast<double>(conform1.size());
    ctr2 /= static_cast<double>(conform2.size());

    Eigen::Matrix3Xd matrix1(3, conform1.size());
    Eigen::Matrix3Xd matrix2(3, conform2.size());
    for (size_t col = 0; col < conform1.size(); col++) {
        matrix1.col(static_cast<int>(col)) = conform1[col] - ctr1;
        matrix2.col(static_cast<int>(col)) = conform2[col] - ctr2;
    }
    matrix2 /= scale;

    // SVD
    Eigen::MatrixXd Cov = matrix1 * matrix2.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Find the rotation
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0) {
        d = 1.0;
    } else {
        d = -1.0;
    }
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;
    Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

    // The final transform
    A.linear() = scale * R;
    A.translation() = ctr2 - scale*R*ctr1;

    return A;
}
