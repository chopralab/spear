// License:

#ifndef SPEAR_SCORINGFUNCTION_HPP
#define SPEAR_SCORINGFUNCTION_HPP

#include <functional>   // std::bad_function_call

#include "spear/exports.hpp"
#include "spear/Molecule.hpp"

namespace Spear {

class SPEAR_EXPORT ScoringFunction {
public:

    virtual ~ScoringFunction() = default;

    /// Calculate the score of a given molecule
    virtual double score(
        const Molecule& mol1,
        const std::vector<size_t>& types1,
        const Molecule& mol2,
        const std::vector<size_t>& types2
    ) = 0;
};

template<class ScoringFunction>
std::string scoringfunction_name() {
    throw std::bad_function_call("scoringfunction_name is unimplemented for this scoring function");
}

}

#endif