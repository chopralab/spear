// License:

#ifndef SPEAR_FORCEFIELD_HPP
#define SPEAR_FORCEFIELD_HPP

#include <unordered_map>
#include <algorithm>

#include <boost/functional/hash.hpp>
#include "spear/exports.hpp"


#include <iostream>
namespace OpenMM {
class System;
}

namespace Spear {

class Molecule;

class SPEAR_EXPORT Forcefield {

public:

    virtual ~Forcefield(){}

    /// Add force to system. `Simulation` will take ownership of the force.
    /// This function must be called after all the particles have been added
    /// to the system.
    virtual void add_forces(const Molecule& mol, OpenMM::System& system) const = 0;

    /// Retrieves the masses of all atoms in a molecule
    virtual std::vector<double> masses(const Molecule& mol) const = 0;

    typedef std::array<size_t, 2> bond_type;
    struct bond_type_hash : public std::unary_function<bond_type, std::size_t> {
        std::size_t operator()(const bond_type& k) const {
            size_t seed = 0;
		    boost::hash_combine(seed, std::min(std::get<0>(k), std::get<1>(k)));
		    boost::hash_combine(seed, std::max(std::get<0>(k), std::get<1>(k)));

            return seed;
        }
    };

    struct bond_type_equal : public std::binary_function<bond_type, bond_type, bool> {
        std::size_t operator()(const bond_type& j, const bond_type& k) const {
            return std::min(k[0], k[1]) == std::min(j[0], j[1]) &&
                   std::max(k[0], k[1]) == std::max(j[0], j[1]);
        }
    };

    typedef std::array<size_t, 3> angle_type;
    struct angle_type_hash : public std::unary_function<angle_type, std::size_t> {
        std::size_t operator()(const angle_type& k) const {
            size_t seed = 0;
		    boost::hash_combine(seed, std::min(std::get<0>(k), std::get<2>(k)));
            boost::hash_combine(seed, std::get<1>(k));
		    boost::hash_combine(seed, std::max(std::get<0>(k), std::get<2>(k)));

            return seed;
        }
    };

    struct angle_type_equal : public std::binary_function<angle_type, angle_type, bool> {
        std::size_t operator()(const angle_type& j, const angle_type& k) const {
            return k[1] == j[1] &&
                   std::min(k[0], k[2]) == std::min(j[0], j[2]) &&
                   std::max(k[0], k[2]) == std::max(j[0], j[2]);
        }
    };

    typedef std::array<size_t, 4> torsion_type;
    struct torsion_type_hash : public std::unary_function<torsion_type, std::size_t> {
        std::size_t operator()(const torsion_type& k) const {
            size_t seed = 0;
            if (std::get<0>(k) < std::get<3>(k)) {
                boost::hash_combine(seed, std::get<0>(k));
                boost::hash_combine(seed, std::get<1>(k));
                boost::hash_combine(seed, std::get<2>(k));
                boost::hash_combine(seed, std::get<3>(k));
            } else if (std::get<0>(k) > std::get<3>(k)) {
                boost::hash_combine(seed, std::get<3>(k));
                boost::hash_combine(seed, std::get<2>(k));
                boost::hash_combine(seed, std::get<1>(k));
                boost::hash_combine(seed, std::get<0>(k));
            } else { // equal
                if (std::get<1>(k) < std::get<2>(k)) {
                    boost::hash_combine(seed, std::get<0>(k));
                    boost::hash_combine(seed, std::get<1>(k));
                    boost::hash_combine(seed, std::get<2>(k));
                    boost::hash_combine(seed, std::get<3>(k));
                } else { // all damn equal or std::get<1>(k) > std::get<2>(k)
                    boost::hash_combine(seed, std::get<3>(k));
                    boost::hash_combine(seed, std::get<2>(k));
                    boost::hash_combine(seed, std::get<1>(k));
                    boost::hash_combine(seed, std::get<0>(k));
                }
            }
            return seed;
        }
    };

    struct torsion_type_equal : public std::binary_function<torsion_type, torsion_type, bool> {
        std::size_t operator()(const torsion_type& j, const torsion_type& k) const {
            if ( (j[0] < j[3] && k[0] < k[3]) || (j[0] > j[3] && k[0] > k[3])) {
                return j == k; // Simple comparison as they are the same order
            }
            if ( (j[0] > j[3] && k[0] < k[3]) || (j[0] < j[3] && k[0] > k[3])) {
                // Reverse order!!!
                return j[0] == k[3] && j[1] == k[2] && j[2] == k[1] && j[3] == k[0];
            }
            if ( j[0] == k[0] && (j[1] < j[2])) {
                // Proper order with ends equal
                return j[1] == std::min(k[1], k[2]) && j[2] == std::max(k[1], k[2]);
            }
            if ( j[0] == k[0] && (j[1] > j[2])) {
                // Reverse order with ends equal
                return j[1] == std::max(k[1], k[2]) && j[2] == std::min(k[1], k[2]);
            }
            // must be all the same
            return j == k;
        }
    };

    typedef std::array<size_t, 4> improper_type;
    struct improper_type_hash : public std::unary_function<improper_type, std::size_t> {
        std::size_t operator()(const improper_type& k) const {
            size_t seed = 0;

            std::array<size_t, 3> others({std::get<0>(k),
                                          std::get<2>(k),
                                          std::get<3>(k)});
            std::sort(others.begin(), others.end());

		    boost::hash_combine(seed, others[0]);
            boost::hash_combine(seed, std::get<1>(k));
		    boost::hash_combine(seed, others[1]);
            boost::hash_combine(seed, others[2]);

            return seed;
        }
    };

    struct improper_type_equal : public std::binary_function<improper_type, improper_type, bool> {
        std::size_t operator()(const improper_type& j, const improper_type& k) const {
            if (j[1] != k[1]) {
                return false;
            }

            std::array<size_t, 3> jother({std::get<0>(j),
                                          std::get<2>(j),
                                          std::get<3>(j)});
            std::sort(jother.begin(), jother.end());

            std::array<size_t, 3> kother({std::get<0>(k),
                                          std::get<2>(k),
                                          std::get<3>(k)});
            std::sort(kother.begin(), kother.end());

            return jother[0] == kother[0] &&
                   jother[1] == kother[1] &&
                   jother[2] == kother[2];
        }
    };
};

}

#endif
