#ifndef SPEAR_CONFIGURATIONS_HPP
#define SPEAR_CONFIGURATIONS_HPP

#include <vector>

#include "Eigen/Geometry"

namespace Spear {

class Configurations {

    using internal_storage = std::list<std::vector<Vector3d>>;

    std::list<std::vector<Vector3d>> configurations_;
    std::vector<Vector3d>* current_configuration_;

public:

    void append_configuration(std::vector<Vector3d> new_config) {
        configurations_
    }

};

}

#endif
