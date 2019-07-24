// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#include "spear/partialcharges/Gasteiger.hpp"
#include "spear/AtomType.hpp"
#include "spear/Graph.hpp"

using namespace Spear;

// Parameters from Table 1 in [Gasteiger 1980]
const std::array<double, 3> parameters[] = {
    {00.00, 00.000, 0.000}, // Skip
    { 7.17,  6.240,-0.560}, // H
    { 7.98,  9.180, 1.880}, // C (sp3)
    {08.79,  9.320, 1.510}, // C (sp2)
    {10.39,  9.450, 0.730}, // C (sp)
    {11.54, 10.820, 1.360}, // N (sp3)
    {00.00, 11.860,11.860}, // N (N3+)
    {12.32, 11.200, 1.340}, // N (N3 planar)
    {12.87, 11.150, 0.850}, // N (sp2)
    {15.68, 11.700,-0.270}, // N (sp)
    {14.18, 12.920, 1.390}, // O (sp3)
    {17.07, 13.790, 0.470}, // O (sp2)
    {14.66, 13.850, 2.310}, // F
    {11.00,  9.690, 1.350}, // Cl
    {10.08,  8.470, 1.160}, // Br
    { 9.90,  7.960, 0.960}, // I
    { 8.90,  8.240, 0.960}, // P (sp3)
    {9.665,  8.530, 0.735}, // P (sp2)
    {10.14,  9.130, 1.380}, // S
    {12.00, 10.805, 1.195}, // S (>SO2)
    { 7.30,  6.567, 0.657}, // Si (sp3)
    {7.905,  6.748, 0.443}, // Si (sp2)
    {9.065,  7.027,-0.002}, // Si (sp)
    { 5.98,  6.820, 1.605}, // B (sp3)
    { 6.42,  6.807, 1.322}, // B (sp2)
    {3.845,  6.755, 3.165}, // Be (sp3)
    {4.005,  6.725, 3.035}, // Be (sp2)
    {3.565,  5.572, 2.197}, // Mg (sp2)
    { 3.30,  5.587, 2.447}, // Mg (sp3)
    { 4.04,  5.472, 1.823}, // Mg (sp)
    {5.375,  4.953, 0.867}, // Al (sp3)
    {5.795,  5.020, 0.695}, // Al (sp2)
};

static size_t get_atom_idx(const AtomVertex& av) {
    switch(av.atomic_number()) {
    case Element::H:
        return 1;
    case Element::C:
        switch (av.hybridization()) {
        case Hybridization::SP3:
            return 2;
        case Hybridization::SP2:
            return 3;
        case Hybridization::SP:
            return 4;
        default:
            return 0;
        }
    case Element::N:
        switch (av.hybridization()) {
        case Hybridization::SP3:
            if (av.degree() >= 4 || av.formal_charge() > 0) {
                return 5;
            }
            if (av.is_planar()) {
                return 6;
            }
            return 7;
        case Hybridization::SP2:
            return 8;
        case Hybridization::SP:
            return 9;
        default:
            return 0;
        }
    case Element::O:
        switch (av.hybridization()) {
        case Hybridization::SP3:
            return 10;
        case Hybridization::SP2:
            return 11;
        default:
            return 0;
        }
    case Element::F:
        return 12;
    case Element::Cl:
        return 13;
    case Element::Br:
        return 14;
    case Element::I:
        return 15;
    case Element::P:
        switch (av.hybridization()) {
        case Hybridization::SP3:
            return 16;
        default:
            return 17;
        }
    case Element::S:
        return 18;
    default:
        return 0;
    }
}

GasteigerCharge::GasteigerCharge(const Molecule& mol) : mol_(mol) {
    resize(mol_.size());
    init_();
}

void GasteigerCharge::init_() {
    std::vector<double> electronegativities(mol_.size());
    std::vector<size_t> atom_index(mol_.size());

    for (auto atom : mol_) {
        auto idx = get_atom_idx(atom);
        atom_index[atom] = idx;
        electronegativities[atom] = parameters[idx][0];
        set(atom, 0.00);
    }

    // run algorithm for six iterations
    for(int iteration = 1; iteration <= 6; iteration++){

        // calculate charges
        for(auto atom : mol_) {
            auto qi = 0.0;
            auto Xi = electronegativities[atom];
            auto pi_index = atom_index[atom];
            if (pi_index == 0) {
                continue;
            }
            const auto& pi = parameters[pi_index];

            for(auto neighbor : atom.neighbors()) {
                auto Xj = electronegativities[neighbor];
                auto pj_index = atom_index[neighbor];
                if (pj_index == 0) {
                    continue;
                }

                const auto& pj = parameters[pj_index];

                auto scale = Xj > Xi ?
                    (atom.atomic_number() == Element::H ? 20.02
                                                        : (pi[0] + pi[1] + pi[2])
                    )
                    :
                    (neighbor.atomic_number() == Element::H ? 20.02
                                                            : (pj[0] + pj[1] + pj[2])
                    );

                scale = 1.0 / scale;
                qi += scale * (Xj - Xi);
            }

            set(atom, get(atom) + qi * pow(0.5, iteration));
        }

        // calculate electronegativities
        for(auto atom : mol_){
            const auto& p = parameters[atom_index[atom]];
            auto Qi = get(atom);

            electronegativities[atom] = p[0] + p[1] * Qi + p[2] * pow(Qi, 2);
        }
    }
}
