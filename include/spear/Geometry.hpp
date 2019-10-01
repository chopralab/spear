// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_GEOMETRY_HPP
#define SPEAR_GEOMETRY_HPP

#include <vector>
#include "spear/exports.hpp"
#include "Eigen/Geometry"

namespace Spear {

using Vector3d = Eigen::Vector3d;
using Conformation = std::vector<Vector3d>;

inline double distance(const Vector3d& u, const Vector3d& v) {
    Vector3d dist = u - v;
    return dist.norm();
}

inline double angle(const Vector3d& a, const Vector3d& b,
                    const Vector3d& c) {
    Vector3d rab = a - b;
    Vector3d rcb = c - b;

    auto cos = rab.dot(rcb) / (rab.norm() * rcb.norm());
    cos = std::max(-1.0, std::min(1.0, cos));
    return std::acos(cos);
}

inline double dihedral(const Vector3d& i, const Vector3d& j,
                       const Vector3d& k, const Vector3d& m) {
    Vector3d rij = i - j;
    Vector3d rjk = j - k;
    Vector3d rkm = k - m;

    Vector3d a = rij.cross(rjk);
    Vector3d b = rjk.cross(rkm);
    return std::atan2(rjk.norm() * b.dot(rij), a.dot(b));
}

inline double nonplanar(const Vector3d& i, const Vector3d& j,
                        const Vector3d& k, const Vector3d& m) {
    Vector3d rji = j - i;
    Vector3d rik = i - k;
    Vector3d rim = i - m;

    Vector3d n = rik.cross(rim);
    return rji.dot(n) / n.norm();
}

inline double square_deviation(const Conformation& c1, const Conformation& c2) {
    if (c1.size() != c2.size()) {
        throw std::runtime_error("Spear::square_deviation(): input data mis-match");
    }

    double distance = 0;
    for (size_t i = 0; i < c1.size(); ++i) {
        distance += (c1[i] - c2[i]).squaredNorm();
    }
    return distance;
}

inline double rmsd(const Conformation& c1, const Conformation& c2) {
    double distance = square_deviation(c1, c2);
    return std::sqrt(distance / static_cast<double>(c1.size()));
}

size_t SPEAR_EXPORT dimensionality(const Conformation& positions, double eps = 0.00001);

Eigen::Affine3d SPEAR_EXPORT kabsch(const Conformation& conform1, const Conformation& conform2);

}

#endif
