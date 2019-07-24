// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_VINASCORE_HPP
#define SPEAR_VINASCORE_HPP

#include "spear/ScoringFunction.hpp"

namespace Spear {

class SPEAR_EXPORT VinaScore final : public ScoringFunction {
public:

    VinaScore(
        double w_g1 = -0.035579,
        double w_g2 = -0.005156,
        double w_rep=  0.840245,
        double w_phobic = -0.035069,
        double w_hydro  = -0.587439,
        double cutoff   = 8.0,
        double w_nrot   = 0.05845
    ) :
        gauss1_weight_(w_g1), gauss2_weight_(w_g2), repulsion_weight_(w_rep),
        hydrophobic_weight_(w_phobic), hydrogen_weight_(w_hydro),
        dist_cutoff_(cutoff), nrotate_weight_(w_nrot) {}

    double score(const Grid& grid, const Molecule& mol1, const Molecule& mol2) const override;
    double score(const Grid& grid, const Molecule& mol, size_t residue_id) const override;
    double score(size_t atomtype1, size_t atomtype2, double r) const override;

    //! An object which represents the components of XScore/Vina's scoring function
    //!
    //! The XScore/Vina scoring function is divided into five components:
    //! two guassian functions, hydrophobic interactions, repulsive forces, and
    //! hydrogen bonding.
    struct Components {
        double g1;          //!< The small gaussian term
        double g2;          //!< The large gaussian term
        double rep;         //!< The repulsive term
        double hydrophobic; //!< The hydrophobic term
        double hydrogen;    //!< The hydrogen bonding term
    };

    //! Calculates the components of a interaction.
    //!
    //! \param atomtype1 The Vina type of an atom
    //! \param atomtype2 The Vina type of an other atom
    //! \param r The distance between the atoms
    Components calculate_components(size_t atomtype1, size_t atomtype2, double r) const;

    //! Calculates the components of a interaction.
    //!
    //! Note: the mol1 and mol2 arguments must have the 'VinaType' atom type
    //! added for this calculation to work, otherwise a Components will be
    //! returned with zero interactions.
    //!
    //! \param grid Grid created from the atom locations in mol1
    //! \param mol1 
    //! \param r The distance between the atoms
    Components calculate_components(const Grid& grid, const Molecule& mol1, const Molecule& mol2) const;

private:
    double gauss1_weight_;
    double gauss2_weight_;
    double repulsion_weight_;
    double hydrophobic_weight_;
    double hydrogen_weight_;
    double dist_cutoff_;
    double nrotate_weight_;
};

template<> std::string SPEAR_EXPORT scoringfunction_name<VinaScore>() {
    return "VinaScore";
}

}

#endif
