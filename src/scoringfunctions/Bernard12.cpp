#include "spear/scoringfunctions/Bernard12.hpp"

#include <fstream>

using namespace Spear;

Bernard12::Bernard12(Options opt, double cutoff,
              const AtomicDistributions& atom_dist,
              const std::unordered_set<size_t>& allowed_atoms):
            options_(opt), dist_cutoff_(cutoff), allowed_atoms_(allowed_atoms) {
    compile_scoring_function_(atom_dist);
}

void Bernard12::compile_scoring_function_(const AtomicDistributions& distributions) {

    step_in_file_ = distributions.step_in_file;
    size_t cutoff_index = static_cast<size_t>(dist_cutoff_ / step_in_file_);

    // Temporary structures to hold the atomic distributions and reductions there of.
    PairVectorDouble gij_of_r; // radial distribution of j given i at distance r
    PairDouble sum_gij_of_r; // sum of above quantity for all distances
    std::vector<double> distance_range_sum(cutoff_index, 0); // sum of all quantities at a distance
    double total_quantity = 0; // sum of all quantities

    // From here on out, pair_dist represents a pair of atoms and their distances
    for (const auto& pair_dist : distributions.values) {
        const auto& atom_pair = pair_dist.first;
        if (options_ & REDUCED &&
            (allowed_atoms_.count(atom_pair.first) == 0 ||
             allowed_atoms_.count(atom_pair.second) == 0)) {
            continue;
        }

        gij_of_r[atom_pair] = std::vector<double>(cutoff_index);
        sum_gij_of_r[atom_pair] = 0;

        for (size_t index = 0; index < pair_dist.second.size(); ++index) {
            auto lower_bound = step_in_file_ * static_cast<double>(index);
            auto upper_bound = step_in_file_ * static_cast<double>(index + 1);

            if (index >= cutoff_index) {
                break;
            }

            auto quantity = pair_dist.second[index];

            auto shell_volume =
                options_ & RADIAL ? 4.0 * M_PI * pow(upper_bound, 3) / 3.0 -
                                    4.0 * M_PI * pow(lower_bound, 3) / 3.0 :
                                    1.0;

            auto shell_density = quantity / shell_volume;

            gij_of_r[atom_pair][index] = shell_density;
            sum_gij_of_r[atom_pair] += shell_density;

            if (options_ & MEAN) continue;
            distance_range_sum[index] += shell_density;
            total_quantity += shell_density;
        }
    }

    if (options_ & MEAN) {
        for (auto& pair_dist : gij_of_r) {

            // Avoid a divide by zero
            if (sum_gij_of_r[pair_dist.first] <= 0) {
                continue;
            }

            for (size_t i = 0; i < pair_dist.second.size(); ++i) {
                distance_range_sum[i] += pair_dist.second[i] /
                                         sum_gij_of_r[pair_dist.first];
            }
        }
    }

    // Now we can actually start compiling the scoring function
    auto num_atom_types = static_cast<double>(options_ & REDUCED ?
        static_cast<size_t>(allowed_atoms_.size()) :
        distributions.max_ids
    );

    // Triangle function
    auto num_pairs = (num_atom_types+1.0)*num_atom_types / 2.0;

    for (auto& pair_dist : gij_of_r) {
        auto w1 = distributions.van_der_waals(pair_dist.first.first);
        auto w2 = distributions.van_der_waals(pair_dist.first.second);
        auto vdW_sum = ((w1 > 0 && w2 > 0) ? w1 + w2 : 4.500);
        auto repulsion_idx = static_cast<size_t>((vdW_sum - 0.6) / step_in_file_);
        auto eps = std::numeric_limits<double>::epsilon();

        energies_scoring_[pair_dist.first] =
            std::vector<double>(pair_dist.second.size(), std::numeric_limits<double>::lowest());

        for (size_t i = 0; i < pair_dist.second.size(); ++i) {

            if (sum_gij_of_r[pair_dist.first] < eps) {
                energies_scoring_[pair_dist.first][i] = 0;
            } else if (pair_dist.second[i] < eps && (i + 1 < repulsion_idx)) {
                energies_scoring_[pair_dist.first][i] = 5.0;
            } else if (pair_dist.second[i] < eps && (i + 1 >= repulsion_idx)) {
                energies_scoring_[pair_dist.first][i] = 0;
            } else {
                auto gij_r_ratio = pair_dist.second[i] /
                                sum_gij_of_r[pair_dist.first];

                auto reference = options_ & MEAN ?
                                 distance_range_sum[i] / num_pairs :
                                 distance_range_sum[i] / total_quantity;

                energies_scoring_[pair_dist.first][i] = -std::log(gij_r_ratio / reference);
            }
        }
    }
}

double Bernard12::score(
    const Molecule& mol1,
    const std::vector<size_t>& types1,
    const Molecule& mol2,
    const std::vector<size_t>& types2
) {

    auto energy_sum = 0.0;
    for (auto atom1 : mol1) {
        if (ignore_hydro && atom1.atomic_number() == 1) continue;
        for (auto atom2 : mol2) {
            if (ignore_hydro && atom2.atomic_number() == 1) continue;
            auto dist = (atom1.position() - atom2.position()).norm();

            auto atom_pair = std::minmax(types1[atom1.index()], types2[atom2.index()]);
            auto idx = static_cast<size_t>(dist / step_in_file_);
            auto distrution = energies_scoring_.find(atom_pair);

            if (distrution == energies_scoring_.end()) {
                continue;
            }

            if (idx < distrution->second.size()) {
                energy_sum += distrution->second[idx];
            }
        }
    }
    return energy_sum;
}
