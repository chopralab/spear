// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include "spear/scoringfunctions/VinaScore.hpp"
#include "spear/atomtypes/VinaType.hpp"
#include "spear/Grid.hpp"
#include "spear/Molecule.hpp"

#include <cmath>

using namespace Spear;

static bool donor_acceptor(size_t t1, size_t t2) {
    return VinaType::is_donor(t1) && VinaType::is_acceptor(t2);
}

static bool h_bond_possible(size_t t1, size_t t2) {
    return donor_acceptor(t1, t2) || donor_acceptor(t2, t1);
}

static double optimal_distance(size_t xs_t1, size_t xs_t2) {
    return van_der_waals<VinaType>(xs_t1) + van_der_waals<VinaType>(xs_t2);
}

static double gaussian(double offset, double width, double r) {
    double exponent = (r - offset) / width;
    return std::exp(-exponent * exponent);
}

static double repulsion(double offset, double r) {
    r -= offset;
    if (r > 0.0) {
        return 0.0;
    }

    return r * r;
}

static double slope_step(double intercept, double end, double r) {
    if(intercept < end) {
        if(r <= intercept) return 0;
        if(r >= end) return 1;
    } else {
        if(r >= intercept) return 0;
        if(r <= end) return 1;
    }
    return (r - intercept) / (end - intercept);
}

double VinaScore::score(const Grid& grid, const Molecule& mol1, const Molecule& mol2) const {
    auto opt_types1 = mol1.atomtype("vina");
    auto opt_types2 = mol2.atomtype("vina");

    if (opt_types1 == nullptr || opt_types2 == nullptr) {
        throw std::invalid_argument("Correct atom types not present in molecule");
    }

    auto& types1 = *opt_types1;
    auto& types2 = *opt_types2;

    Components comps = calculate_components(grid, mol1, mol2);
    return comps.g1          * gauss1_weight_ +
           comps.g2          * gauss2_weight_ +
           comps.rep         * repulsion_weight_ +
           comps.hydrophobic * hydrophobic_weight_ +
           comps.hydrogen    * hydrogen_weight_;
}

double VinaScore::score( const Grid& grid, const Molecule& mol, size_t residue_id) const {
    auto opt_types1 = mol.atomtype("vina");
    auto& types1 = *opt_types1;

    const auto& residue = mol.topology().residues()[residue_id];

    auto energy_sum = 0.0;
    for (auto atom2_id : residue) {
        auto atom2 = mol[atom2_id];
        if (atom2.name() == "C" || atom2.name() == "N" || atom2.name() == "O" || atom2.name() =="CA") continue;
        auto neighbors = grid.neighbors(atom2.position(), dist_cutoff_);
        for (auto neighbor : neighbors) {
            auto atom1 = mol[neighbor];
            if (residue.contains(neighbor)) continue;
            auto dist = (atom1.position() - atom2.position()).norm();

            energy_sum += score(types1[atom1], types1[atom2], dist);
        }
    }
    return energy_sum;
}

double VinaScore::score(size_t atomtype1, size_t atomtype2, double r) const {
    if (r > dist_cutoff_) {
        return 0.0;
    }

    auto comps = calculate_components(atomtype1, atomtype2, r);
    return comps.g1          * gauss1_weight_ +
           comps.g2          * gauss2_weight_ +
           comps.rep         * repulsion_weight_ +
           comps.hydrophobic * hydrophobic_weight_ +
           comps.hydrogen    * hydrogen_weight_;
}

VinaScore::Components VinaScore::calculate_components(size_t atomtype1, size_t atomtype2, double r) const {
    VinaScore::Components ret{0.0, 0.0, 0.0, 0.0, 0.0};

    if (VinaType::is_skip(atomtype1) ||
        VinaType::is_skip(atomtype2)) {
        return ret;
    }

    r -= optimal_distance(atomtype1, atomtype2);
    ret.g1 = gaussian(0.0, 0.5, r);
    ret.g2 = gaussian(3.0, 2.0, r);

    ret.rep = repulsion(0, r);
    if (VinaType::is_hydrophobic(atomtype1) &&
        VinaType::is_hydrophobic(atomtype2)) {
        ret.hydrophobic = slope_step(1.5, 0.5, r);
    }

    if (h_bond_possible(atomtype1, atomtype2)) {
        ret.hydrogen = slope_step(0.0, -0.7, r);
    }

    return ret;
}

VinaScore::Components VinaScore::calculate_components(const Grid& grid, const Molecule& mol1, const Molecule& mol2) const {
    auto opt_types1 = mol1.atomtype("vina");
    auto opt_types2 = mol2.atomtype("vina");

    if (opt_types1 == nullptr || opt_types2 == nullptr) {
        throw std::invalid_argument("Correct atom types not present in molecule");
    }

    auto& types1 = *opt_types1;
    auto& types2 = *opt_types2;

    VinaScore::Components ret{0.0, 0.0, 0.0, 0.0, 0.0};
    for (auto atom2 : mol2) {
        auto neighbors = grid.neighbors(atom2.position(), dist_cutoff_);
        for (auto neighbor : neighbors) {
            auto atom1 = mol1[neighbor];
            auto dist = (atom1.position() - atom2.position()).norm();
            if (dist > dist_cutoff_) {
                continue;
            }

            auto component = calculate_components(types1[atom1], types2[atom2], dist);
            ret.g1 += component.g1;
            ret.g2 += component.g2;
            ret.rep += component.rep;
            ret.hydrogen += component.hydrogen;
            ret.hydrophobic += component.hydrophobic;
        }
    }

    return ret;
}
