#ifndef SPEAR_BERNARD12_HPP
#define SPEAR_BERNARD12_HPP

#include "spear/ScoringFunction.hpp"
#include "spear/AtomicDistributions.hpp"

namespace Spear {

class SPEAR_EXPORT Bernard12 : public ScoringFunction {
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
              const AtomicDistributions& atom_dist, const std::string& atomtype,
              const std::unordered_set<size_t>& allowed_atoms = std::unordered_set<size_t>());

    double score(const Molecule& mol1, const Molecule& mol2) override;

    bool ignore_hydro = false;
private:

    PairVectorDouble energies_scoring_;  // scoring function

    Options options_;
    double dist_cutoff_;
    double step_in_file_;
    std::unordered_set<size_t> allowed_atoms_;
    const std::string atomtype_;

    void compile_scoring_function_(const AtomicDistributions& distributions);
};

template<> std::string scoringfunction_name<Bernard12>();

}

#endif
