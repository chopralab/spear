// License:

#ifndef SPEAR_SCORINGFUNCTION_HPP
#define SPEAR_SCORINGFUNCTION_HPP

#include <functional>   // std::bad_function_call
#include <string>

#include "spear/exports.hpp"

namespace Spear {

class Molecule;
class Grid;

class SPEAR_EXPORT ScoringFunction {
public:

    virtual ~ScoringFunction() = default;

    /// Calculate the score of a given molecule
    virtual double score( const Grid& grid, const Molecule& mol1, const Molecule& mol2) const = 0;

    virtual double score( const Grid& grid, const Molecule& mol, size_t residue_id) const = 0;
};

template<class ScoringFunction>
std::string scoringfunction_name() {
    throw std::bad_function_call();
}

}

#endif
