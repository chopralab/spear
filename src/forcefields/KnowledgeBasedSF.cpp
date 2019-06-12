#include "spear/forcefields/KnowledgeBasedSF.hpp"
#include "spear/Molecule.hpp"
#include "spear/atomtypes/IDATM.hpp"
#include "OpenMM.h"

#include <unordered_set>
#include <numeric>

using namespace Spear;

auto IDATM_name = atomtype_id_for_name<IDATM>;

static std::pair<size_t, size_t> name_pair(const std::string& n1, const std::string& n2) {
    auto t1 = IDATM_name(n1);
    auto t2 = IDATM_name(n2);
    auto mm = std::minmax(t1, t2);
    return std::make_pair(mm.first, mm.second);
}

void KnowledgeBasedSF::add_forces(const std::vector<std::reference_wrapper<const Molecule>>& mols, OpenMM::System& system) const {
    using namespace OpenMM;

    auto forcefield = new CustomNonbondedForce(
        "scale * ( step(rep_idx - r) * (B*(0.1/r)^12) + kbpot(r, idatm1 , idatm2) ); B = btable(idatm1, idatm2); rep_idx = rtable(idatm1, idatm2);"
    );

    forcefield->setNonbondedMethod(
        CustomNonbondedForce::CutoffNonPeriodic
    );

    forcefield->setCutoffDistance(1.5);
    forcefield->addGlobalParameter("scale", 1);
    forcefield->addPerParticleParameter("idatm");

    std::vector<std::pair<int, int>> bond_pairs;
    std::map<size_t, size_t> idatm_to_internal;
    std::vector<size_t> used_atom_types;
    size_t num_types = 0;
    size_t num_atoms = 0;

    for (auto mol : mols) {
        auto at = mol.get().atomtype(atomtype_);
        if (at == nullptr) {
            //TODO: throw something
        }

        for (auto atom : mol.get()) {
            auto atype = (*at)[atom];
            if (!idatm_to_internal.count(atype)) {
                idatm_to_internal[atype] = num_types;
                used_atom_types.push_back(atype);
                num_types++;
            }

            forcefield->addParticle(
                {static_cast<double>(idatm_to_internal[atype])}
            );
        }

        for (auto& bond : mol.get().topology().bonds()) {
            bond_pairs.emplace_back(
                bond[0] + num_atoms,
                bond[1] + num_atoms
            );
        }

        num_atoms += mol.get().size();
    }

    forcefield->createExclusionsFromBonds(bond_pairs, 4);

    std::vector<double> all_energies;
    all_energies.reserve(num_types * num_types * 150);

    std::vector<double> btable;
    btable.reserve(num_types * num_types);

    std::vector<double> rtable;
    rtable.reserve(num_types * num_types);

    for (auto at1 : used_atom_types) {
        for (auto at2 : used_atom_types) {
            double B, rep_idx;
            auto energies = obtain_scores_for_pair_(at1, at2, B, rep_idx);

            all_energies.insert(
                std::end(all_energies),
                std::begin(energies), std::end(energies)
            );

            btable.push_back(B);
            rtable.push_back(rep_idx * 0.1);
        }
    }

    forcefield->addTabulatedFunction(
        "kbpot",
        new Continuous3DFunction(150, num_types, num_types,
                                 all_energies,
                                 0.0, (150 - 1) * 0.01,
                                 0.0, num_types - 1.0,
                                 0.0, num_types - 1.0
        )
    );

    forcefield->addTabulatedFunction(
        "btable",
        new Discrete2DFunction(num_types, num_types,
                               btable
        )
    );

    forcefield->addTabulatedFunction(
        "rtable",
        new Discrete2DFunction(num_types, num_types,
                               rtable
        )
    );

    system.addForce(forcefield);
}

std::vector<double> KnowledgeBasedSF::obtain_scores_for_pair_(size_t at1, size_t at2, double& B, double& rep_index) const {
    
    const std::map<std::pair<size_t, size_t>, double> rep_indicies {
        {name_pair("F",   "Sar"), 3.6},
        {name_pair("Npl", "Npl"), 2.5},
        {name_pair("Npl", "S3"), 2.8},
        {name_pair("O2",  "O2-"), 2.2},
        {name_pair("O2-", "S3"), 2.9},
        {name_pair("Cac", "S3"), 3.4},
        {name_pair("Cl",  "O2-"), 2.8},
        {name_pair("N2",  "O2-"), 2.8},
        {name_pair("N2",  "S3"), 3.0},
        {name_pair("N3+", "N3+"), 2.5},
        {name_pair("N3+", "Ng+"), 2.9},
        {name_pair("N3+", "O2-"), 2.5},
        {name_pair("N3+", "O3"), 2.5},
        {name_pair("N3+", "S3"), 3.0},
        {name_pair("Ng+", "Ng+"), 2.7}
    };

    auto vdw_sum = van_der_waals<IDATM>(at1) + van_der_waals<IDATM>(at1);

    std::vector<double> energy;
    energy.reserve(150);

    for (auto r = 0.0; r + 0.0001 < 15.0; r += 0.1) {
        auto energy_at_point = sf_.score(at1, at2, r + 0.0001);
        energy.push_back(energy_at_point);
    }

    // Find values in potential returned by the scoring function and correct
    // potential outliers
    auto start_idx = static_cast<size_t>((vdw_sum - 0.6) / 0.1);
    for (auto i = start_idx; i < energy.size() - 2; ++i) {
        if (energy[i] == 0.0 || energy[i] == 5.0) { // Skip the 'fake' potential
            continue;
        }

        auto j = i + 1;
        while (j < i + 3 && j < energy.size() &&
                (energy[j] == 0.0 || energy[j] == 5.0)) {
            ++j;
        }
        if (j <= i + 1 || j >= energy.size()) {
            continue;
        }
        auto k0 = (energy[j] - energy[i]) / static_cast<double>(j - i);

        // correct points between 
        for (size_t k = 1; k < j - i; ++k) {
            energy[i + k] = energy[i] + k0 * k;
        }
    }

    // identify the repulsive part of the potential.
    size_t repulsion_idx = 0;
    auto repulsive_look_up = rep_indicies.find(std::minmax(at1, at2));

    if (repulsive_look_up != rep_indicies.end()) {
        repulsion_idx = static_cast<size_t>(repulsive_look_up->second / 0.1);
    } else {
        double global_min = 10000000000.0;
        auto end_idx = static_cast<size_t>((vdw_sum + 1.0) / 0.1);
        for (auto i = start_idx; i < end_idx && i < energy.size(); ++i) {
            if (energy[i] >= global_min) {
                continue;
            }

            global_min = energy[i];
            repulsion_idx = i;

            // ensure the derivative is negative and large 
            while( repulsion_idx >= 1 && energy[repulsion_idx] != 0.0 &&
                ((energy[repulsion_idx] - energy[repulsion_idx - 1]) / 0.1) < 0.75) {
                --repulsion_idx;
            }
        }
    }

    // Calculate 4 derivatives after the repulsive potential
    std::vector<double> derivative;
    for (auto i = repulsion_idx; i < repulsion_idx + 5 && i < energy.size() - 1; ++i) {
        derivative.push_back((energy[i + 1] - energy[i]) / 0.1);
    }

    std::set<size_t> lowest_derivate_indicies;
	for (size_t i = 0; i < 3 && i < derivative.size(); ++i) {
		auto it = std::min_element(derivative.begin(), derivative.end());
        if (*it >= 0.0) { // Don't allow positive derivatives
            break;
        }
		auto i0 = static_cast<size_t>(repulsion_idx + (it - derivative.begin()));
		derivative.erase(it);
		lowest_derivate_indicies.insert(i0);
		lowest_derivate_indicies.insert(i0 + 1);
	}

    std::vector<double> r, p;
    auto p_min = 100000.0;
    if (lowest_derivate_indicies.size() < 2) {
        auto sum = std::accumulate(energy.begin(), energy.end(), 0);
        if (sum == 0) {
            // It's all zeros anyway
            return energy;
        }
        std::cerr << "The potential does not have a repulsion index! Found: "
                  << repulsion_idx << "(" << energy[repulsion_idx] << ")" << std::endl;
        //while (repulsion_idx < energy.size() && energy[repulsion_idx] != 0.0) {
        //    repulsion_idx++;
        //}

        return energy;
    } else {
        repulsion_idx = *(lowest_derivate_indicies.begin());
        auto rep_end = *(lowest_derivate_indicies.rbegin());

        for (auto i = repulsion_idx; i <= rep_end; ++i) {
            r.push_back(i * 0.1);
            p.push_back(energy[i]);
            p_min = std::min(energy[i], p_min);
        }
    }

    auto xsum = 0.0;
    for (auto xs : r) {
        xsum += std::log(xs);
    }

    auto ysum = 0.0;
    for (auto ys : p) {
        ysum += std::log(ys - p_min + 1.00); // Ensure it will be > 0.0
    }

    B = std::exp((ysum + 12.0 * xsum) / static_cast<double>(r.size()));

    auto cost = 0.0;
    for (size_t j = 0; j < r.size(); ++j) {
        cost += B / std::pow(r[j], 12) - p[j];
    }
    cost = cost /  static_cast<double>(r.size());
/*
    std::vector<double> repulsion;
    repulsion.reserve(repulsion_idx);
    auto x1 = repulsion_idx * 0.1;
    for (double xi = 0.1; xi < x1; xi += 0.1) {
        repulsion.push_back(B / std::pow(xi, 12) - cost);
    }
*/
    for (size_t j = 0; j < repulsion_idx; ++j) {
        energy[j] = 0;
    }

    rep_index = static_cast<double>(repulsion_idx) * 0.1;

    return energy;
}
