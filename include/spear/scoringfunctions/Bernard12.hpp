// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_BERNARD12_HPP
#define SPEAR_BERNARD12_HPP

#include "spear/ScoringFunction.hpp"

#include <unordered_set>

#include "spear/AtomicDistributions.hpp"

namespace Spear {

struct AtomicDistributions;

class SPEAR_EXPORT Bernard12 final : public ScoringFunction {
public:

    enum Options {
        RADIAL = 1,
        NORMALIZED_FREQUENCY = 2,
        MEAN = 4,
        CUMULATIVE = 8,
        REDUCED = 16,
        COMPLETE = 32,
    };

    Bernard12(Options opt, double cutoff,
              const AtomicDistributions& atom_dist, std::string atomtype,
              std::unordered_set<size_t> allowed_atoms = std::unordered_set<size_t>());

    double score(const Grid& grid, const Molecule& mol1, const Molecule& mol2) const override;
    double score(const Grid& grid, const Molecule& mol, size_t residue_id) const override;
    double score( size_t atomtype1, size_t atomtype2, double r) const override;

    bool ignore_hydro = true;
private:

    PairVectorDouble energies_scoring_;  // scoring function

    Options options_;
    double dist_cutoff_;
    double step_in_file_;
    std::unordered_set<size_t> allowed_atoms_;
    const std::string atomtype_;

    void compile_scoring_function_(const AtomicDistributions& distributions);
};

template<> std::string SPEAR_EXPORT scoringfunction_name<Bernard12>();

}

#endif
