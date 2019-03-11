#ifndef SPEAR_GEOMETRY_HPP
#define SPEAR_GEOMETRY_HPP

#include "Eigen/Geometry"

namespace Spear {

using Vector3d = Eigen::Vector3d;

double distance(const Vector3d& u, const Vector3d& v) {
    Vector3d dist = u - v;
    return dist.norm();
}

double angle(const Vector3d& a, const Vector3d& b,
             const Vector3d& c) {
    Vector3d rab = a - b;
    Vector3d rcb = c - b;

    auto cos = rab.dot(rcb) / (rab.norm() * rcb.norm());
    cos = std::max(-1.0, std::min(1.0, cos));
    return std::acos(cos);
}

double dihedral(const Vector3d& i, const Vector3d& j,
                const Vector3d& k, const Vector3d& m) {
    Vector3d rij = i - j;
    Vector3d rjk = j - k;
    Vector3d rkm = k - m;

    Vector3d a = rij.cross(rjk);
    Vector3d b = rjk.cross(rkm);
    return atan2(rjk.norm() * b.dot(rij), a.dot(b));
}

double improper(const Vector3d& i, const Vector3d& j,
                const Vector3d& k, const Vector3d& m) {
    Vector3d rji = j - i;
    Vector3d rik = i - k;
    Vector3d rim = i - m;

    Vector3d n = rik.cross(rim);
    return rji.dot(n) / n.norm();
}

}

#endif
