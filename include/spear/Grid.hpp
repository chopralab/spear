// License:

#ifndef SPEAR_GRID_HPP
#define SPEAR_GRID_HPP

#include "spear/exports.hpp"

#include <unordered_map>
#include <unordered_set>

#include "spear/Molecule.hpp"

#include <boost/functional/hash.hpp>

namespace Spear {

using Eigen::Vector3d;

struct SPEAR_EXPORT GridPoint : public std::array<uint64_t, 3> {

    GridPoint(const Vector3d& point,
              const Vector3d& min,
              const Vector3d& max,
              const Vector3d& step,
              double dist = 0.0) {

        // Do not use auto, it causes a strange bug here as it attempts to
        // modify const memory?
        Vector3d shifted = point - min;
        shifted += Vector3d(dist, dist, dist);

        auto xcrd = shifted[0] < 0.0    ? 0.0 :
                    shifted[0] > max[0] ? max[0] :
                    std::floor(shifted[0] / step[0]);

        auto ycrd = shifted[1] < 0.0    ? 0.0 :
                    shifted[1] > max[1] ? max[1] :
                    std::floor(shifted[1] / step[1]);

        auto zcrd = shifted[2] < 0.0    ? 0.0 :
                    shifted[2] > max[2] ? max[2] :
                    std::floor(shifted[2] / step[2]);

        (*this)[0] = static_cast<uint64_t>(xcrd);
        (*this)[1] = static_cast<uint64_t>(ycrd);
        (*this)[2] = static_cast<uint64_t>(zcrd);
    }

	friend size_t hash_value(const GridPoint& p) {
		size_t seed = 0;
		boost::hash_combine(seed, p[0]);
		boost::hash_combine(seed, p[1]);
		boost::hash_combine(seed, p[2]);

		return seed;
	}

    GridPoint(uint64_t i, uint64_t j, uint64_t k)
        : std::array<uint64_t, 3>({i, j, k}) {}
};

inline bool operator==(const GridPoint& a, const GridPoint& b) {
	return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

class SPEAR_EXPORT Grid {
    constexpr static auto eps = std::numeric_limits<double>::epsilon();
    constexpr static auto max = std::numeric_limits<double>::max();
    constexpr static auto min = std::numeric_limits<double>::min();
public:
    Grid(const std::vector<Vector3d>& points, const Vector3d& step = {1.0, 1.0, 1.0})
        : min_(max,max,max), max_(min,min,min), step_(step) {

        if (step_[0] <= eps || step_[1] <= eps || step_[2] <= eps) {
            throw std::invalid_argument("Step size too small");
        }

        for (auto point : points) {
            min_[0] = std::min(point[0], min_[0]);
            min_[1] = std::min(point[1], min_[1]);
            min_[2] = std::min(point[2], min_[2]);

            max_[0] = std::max(point[0], max_[0]);
            max_[1] = std::max(point[1], max_[1]);
            max_[2] = std::max(point[2], max_[2]);
        }

        // Compute the maximum range
        max_ -= min_;

        for (size_t i = 0; i < points.size(); ++i) {
            grid_.emplace(GridPoint(points[i], min_, max_, step_), i);
        }
    }

    std::unordered_set<size_t> neighbors(const Vector3d& point) const {
        auto gridpoint = GridPoint(point, min_, max_, step_);
        std::unordered_set<size_t> ret;

        auto range = grid_.equal_range(gridpoint);
        for (auto gp = range.first; gp != range.second; ++gp) {
            ret.insert(gp->second);
        }

        return ret;
    }

    std::unordered_set<size_t> neighbors(const Vector3d& point,
                                         double dist,
                                         double lower_tol = 0.0,
                                         double upper_tol = 0.0) const {
        auto cmin = GridPoint(point, min_, max_, step_, -(dist + lower_tol));
        auto cmax = GridPoint(point, min_, max_, step_,   dist + upper_tol);
        std::unordered_set<size_t> ret;

        for (auto i = cmin[0]; i <= cmax[0]; ++i) {
            for (auto j = cmin[1]; j <= cmax[1]; ++j) {
                for (auto k = cmin[2]; k <= cmax[2]; ++k) {
                    auto range = grid_.equal_range(GridPoint(i, j, k));
                    for (auto gp = range.first; gp != range.second; ++gp) {
                        ret.insert(gp->second);
                    }
                }
            }
        }

        return ret;
    }

    size_t size() const {
        return grid_.size();
    }

    std::unordered_set<GridPoint, boost::hash<GridPoint>> occupied() const {
        std::unordered_set<GridPoint, boost::hash<GridPoint>> ret;

        for (auto gp : grid_) {
            ret.insert(gp.first);
        }

        return ret;
    }

private:
    Vector3d min_;
    Vector3d max_;
    Vector3d step_;

    std::unordered_multimap<GridPoint, size_t, boost::hash<GridPoint>> grid_;
};

}

#endif
