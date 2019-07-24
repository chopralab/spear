// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_ATOMICDISTRIBUTIONS_HPP
#define SPEAR_ATOMICDISTRIBUTIONS_HPP

#include <map>
#include <vector>
#include <fstream>
#include <algorithm>

#include "spear/AtomType.hpp"

namespace Spear {

typedef std::map<std::pair<size_t, size_t>, std::vector<double>> PairVectorDouble;
typedef std::map<std::pair<size_t, size_t>, double> PairDouble;

struct AtomicDistributions {
    PairVectorDouble values;
    size_t max_ids;
    double step_in_file;
    double max_distance;
    std::function<double(size_t)> van_der_waals;
};

template <class Type>
AtomicDistributions read_atomic_distributions(std::ifstream& atomic_dist) {

    AtomicDistributions distribuitons;
    distribuitons.max_ids = atomtype_id_count<Type>();
    distribuitons.van_der_waals = van_der_waals<Type>;

    atomic_dist >> distribuitons.step_in_file >> distribuitons.max_distance;
    auto number_of_lines_per_pair = static_cast<size_t>(
        distribuitons.max_distance / distribuitons.step_in_file
    );

    while (!atomic_dist.eof()) {

        std::string idatm1, idatm2, tag;
        atomic_dist >> idatm1 >> idatm2 >> tag;

        if (idatm1 == "" || idatm2 == "") {
            break;
        }

        auto idatm1_id = atomtype_id_for_name<Type>(idatm1);
        auto idatm2_id = atomtype_id_for_name<Type>(idatm2);
        auto atom_pair = std::minmax(idatm1_id, idatm2_id);

        if (tag == "null") {
            continue;
        }

        for (size_t i = 0; i < number_of_lines_per_pair; ++i) {
            double count;
            atomic_dist >> count;

            auto map_loc = distribuitons.values.find(atom_pair);
            if (map_loc == distribuitons.values.end()) {
                map_loc = distribuitons.values.emplace(
                    atom_pair,
                    std::vector<double>(number_of_lines_per_pair)
                ).first;
            }

            map_loc->second[i] = count;
        }
    }

    return distribuitons;
}

}

#endif
