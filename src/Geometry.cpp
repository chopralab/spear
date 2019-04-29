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
