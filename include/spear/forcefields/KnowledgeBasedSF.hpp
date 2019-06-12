// License:

#ifndef SPEAR_KNOWNLEDGEBASEDSF_FORCEFIELD_HPP
#define SPEAR_KNOWNLEDGEBASEDSF_FORCEFIELD_HPP

#include "spear/Forcefield.hpp"
#include "spear/ScoringFunction.hpp"

namespace Spear {

class AtomVertex;

class SPEAR_EXPORT KnowledgeBasedSF : public NonBondedForcefield {
public:
    KnowledgeBasedSF(const ScoringFunction& sf, const std::string& atomtype): sf_(sf), atomtype_(atomtype) {}

    void add_forces(const std::vector<std::reference_wrapper<const Molecule>>& mols, OpenMM::System& system) const override;

//private:
    std::vector<double> obtain_scores_for_pair_(size_t at1, size_t at2, double& B, double& rep_index) const;
    const ScoringFunction& sf_;
    std::string atomtype_;
};

}

#endif
