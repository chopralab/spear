// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include <array>
#include <memory>
#include <vector>
#include <algorithm>
#include "spear/Geometry.hpp"

namespace Spear {

struct Bounds
{
    Vector3d center;
    double radius;
};

// ----------------------------------------------------------------------------
// Based on http://www.flipcode.com/archives/Octree_Implementation.shtml
// ----------------------------------------------------------------------------

class Octree
{
public:
    Octree(const Conformation& conf, size_t threshold, size_t max_depth);
    Octree();

    virtual ~Octree() = default;

    static Bounds bounds(const Conformation& points);

    size_t size() const;

    std::vector<size_t> neighbors(const Vector3d& point);

protected:
    void build(
        std::vector<size_t>::const_iterator begin,
        std::vector<size_t>::const_iterator end,
        const Conformation& conf,
        size_t threshold,
        size_t maximum_depth,
        const Bounds &bounds,
        size_t current_depth = 0
    );

    std::array<std::unique_ptr<Octree>, 8> child_;
    std::vector<size_t> points_;
    Vector3d center_;
    double radius_;
};

Octree::Octree(const Conformation& conf, size_t threshold, size_t max_depth) :Octree() {
    std::vector<size_t> points;
    points.reserve(conf.size());
    for (size_t i = 0; i < conf.size(); ++i) {
        points.push_back(i);
    }
    auto bounds = Octree::bounds(conf);
    build(
        std::begin(points), std::end(points),
        conf, threshold, max_depth, bounds, 0
    );
}

Octree::Octree() : center_{0,0,0}, radius_(0.0),
    child_ {nullptr, nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr}
{
}

inline void Octree::build(
    std::vector<size_t>::const_iterator begin,
    std::vector<size_t>::const_iterator end,
    const Conformation& conf,
    size_t threshold,
    size_t maximum_depth,
    const Bounds &bounds,
    size_t current_depth) {

    // Are we a leaf?
    auto count = std::distance(begin, end);
    if (count <= threshold || current_depth >= maximum_depth) {
        // This is a leaf, copy the points over
        std::copy(begin, end, std::back_inserter(points_));
        return;
    }

    std::array<size_t, 8> child_point_counts = {0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<size_t> codes;
    codes.reserve(count);

    center_ = bounds.center;

    // Classify each point to a child node (IE spacial partition)
    std::for_each(begin, end, [&](size_t p){

        std::size_t code = 0;
        if (conf[p][0] > center_[0]) code |= 1;
        if (conf[p][1] > center_[1]) code |= 2;
        if (conf[p][2] > center_[2]) code |= 4;

        child_point_counts[code]++;
        codes.emplace_back(code);
    });

    // Create children for each spatial partition
    for (int i = 0; i < 8; i++) {
        // Blank spatial partitions are nullptr
        if (!child_point_counts[i]) continue;

        child_[i] = std::make_unique<Octree>();

        std::vector<size_t> newList;
        newList.reserve(child_point_counts[i]);
        std::size_t current_copy = 0;
        for (auto code : codes) {
            if (code == i) {
                newList.push_back(*(begin + current_copy));
            }
            current_copy++;
        }

        const Vector3d offset_table[8] =
        {
            {-0.5, -0.5, -0.5},
            {+0.5, -0.5, -0.5},
            {-0.5, +0.5, -0.5},
            {+0.5, +0.5, -0.5},
            {-0.5, -0.5, +0.5},
            {+0.5, -0.5, +0.5},
            {-0.5, +0.5, +0.5},
            {+0.5, +0.5, +0.5}
        };

        auto offset = offset_table[i] * bounds.radius;

        // Each sub space partition has half the radius of the current partition
        // and has its center offset given from the table above
        Bounds new_bounds;
        new_bounds.radius = bounds.radius * 0.5;
        new_bounds.center = bounds.center + offset;

        child_[i]->build(
            std::begin(newList),
            std::end(newList),
            conf,
            threshold, maximum_depth,
            new_bounds, current_depth + 1
        );
    }
}

inline Bounds Octree::bounds(const Conformation& points) {

    auto min = points[0];
    auto max = points[0];

    for (const auto& p : points) {
        if (p[0] < min[0]) min[0] = p[0];
        if (p[1] < min[1]) min[1] = p[1];
        if (p[2] < min[2]) min[2] = p[2];
        if (p[0] > max[0]) max[0] = p[0];
        if (p[1] > max[1]) max[1] = p[1];
        if (p[2] > max[2]) max[2] = p[2];
    }

    // The radius in each direction and its use to find the center
    auto radius = max - min;

    Bounds res;
    res.center = min + radius * 0.5;

    // Use the largest radius
    res.radius = radius[0];
    if (res.radius < radius[1]) res.radius = radius[1];
    if (res.radius < radius[2]) res.radius = radius[2];

    return res;
}

inline size_t Octree::size() const {
    auto size = points_.size();

    if (size != 0) {
        return size;
    }

    for (size_t i = 0; i < 8; ++i) {
        if (child_[i]) {
            size += child_[i]->size();
        }
    }
    return size;
}

inline std::vector<size_t> Octree::neighbors(const Vector3d& point) {
    if (points_.size()) {
        return points_;
    }

    std::size_t code = 0;
    if (point[0] > center_[0]) code |= 1;
    if (point[1] > center_[1]) code |= 2;
    if (point[2] > center_[2]) code |= 4;

    if (!child_[code]) {
        return std::vector<size_t>();
    }

    return child_[code]->neighbors(point);
}

}
