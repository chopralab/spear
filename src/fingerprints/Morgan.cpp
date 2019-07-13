#include "spear/Molecule.hpp"
#include "spear/fingerprints/Morgan.hpp"

#include <numeric>
#include <boost/dynamic_bitset.hpp>

#include <iostream>

using namespace Spear;

using BitSet = std::vector<bool>;
using AccumTuple = std::tuple<boost::dynamic_bitset<>, size_t, unsigned int>;
using BitInfoMap = std::map<std::uint32_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>>;

namespace Spear {

static std::vector<size_t> atom_invarients(const Molecule& mol, bool use_rings) {
    std::vector<size_t> invars;
    invars.reserve(mol.size());

    for (auto atom : mol) {
        size_t seed = 0;
        boost::hash_combine(seed, static_cast<uint64_t>(atom.atomic_number()));
        boost::hash_combine(seed, atom.expected_bonds());
        boost::hash_combine(seed, atom.total_hydrogens());
        boost::hash_combine(seed, atom.formal_charge());
        /* int deltaMass = static_cast<int>(
            atom.mass() -
            PeriodicTable::getTable()->geatomicWeight(atom->geatomicNum()));
        components.push_back(deltaMass);*/
        boost::hash_combine(seed, static_cast<uint32_t>(0));

        const auto& sssrs = atom.sssrs();

        if (use_rings && sssrs.first != sssrs.second) {
            boost::hash_combine(seed, 1ull);
        }
        invars.push_back(seed);
    }

    return invars;
}

template<typename Container>
static bool contains(const Container& c, const typename Container::value_type& m) {
    return std::find(std::begin(c), std::end(c), m) != std::end(c);
}

Fingerprint
calcFingerprint(const Molecule& mol, uint64_t radius,
                const std::vector<size_t>& invariants_in,
                const std::vector<size_t>& fromAtoms,
                size_t finger_print_length,
                bool useChirality, bool useBondTypes,
                bool onlyNonzeroInvariants) {

    auto nAtoms = mol.size();

    BitInfoMap atomsSettingBits;
    Fingerprint res(finger_print_length);

    std::vector<size_t> invariants;

    if (invariants_in.empty()) {
        invariants = atom_invarients(mol, true);
    } else {
        invariants = invariants_in;
    }

    // add the round 0 invariants to the result:
    for (size_t i = 0; i < nAtoms; ++i) {
        if (!fromAtoms.empty() && !contains(fromAtoms, i)) {
            continue;
        }

        if (onlyNonzeroInvariants && invariants[i] == 0) {
            continue;
        }

        auto bit = res.increase(invariants[i]);
        atomsSettingBits[bit].push_back(std::make_pair(i, 0));
    }

    BitSet chiralAtoms(nAtoms);

    // these are the neighborhoods that have already been added
    // to the fingerprint
    std::vector<boost::dynamic_bitset<>> neighborhoods;

    // these are the environments around each atom:
    auto atom_neighborhood = std::vector<boost::dynamic_bitset<>>(
        nAtoms, boost::dynamic_bitset<>(mol.bond_count())
    );

    // atoms with exact environments already seen
    BitSet deadAtoms(nAtoms, false);

    BitSet include_atom(nAtoms, false);
    if (!fromAtoms.empty()) {
        for (auto idx : fromAtoms) {
            include_atom[idx] = true;
        }
    } else {
        include_atom = BitSet(nAtoms, true);
    }

    std::vector<size_t> atom_order;
    if (onlyNonzeroInvariants) {
        std::vector<std::pair<size_t, size_t>> ordering(nAtoms);

        std::generate(ordering.begin(), ordering.end(),
            [&invariants, i = 0] () mutable {
                if (invariants[i] == 0) {
                    return std::make_pair(1, i);
                }
                return std::make_pair(0, i);
            }
        );

        std::sort(ordering.begin(), ordering.end());
        atom_order.reserve(nAtoms);
        std::transform(ordering.cbegin(), ordering.cend(),
            std::back_inserter(atom_order),
            [](const std::pair<size_t, size_t>& other){
                return other.second;
            }
        );
    } else {
        atom_order.resize(nAtoms);
        std::iota(atom_order.begin(), atom_order.end(), 0);
    }

    // now do our subsequent rounds:
    for (uint64_t layer = 0; layer < radius; ++layer) {
        std::vector<size_t> roundInvariants(nAtoms);
        std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods = atom_neighborhood;
        std::vector<AccumTuple> round_neighborhoods;

        for (auto atomIdx : atom_order) {
            if (deadAtoms[atomIdx]) {
                continue;
            }

            auto atom = mol[atomIdx];
            if (atom.degree() == 0) {
                deadAtoms[atomIdx] = true;
                continue;
            }

            std::vector<std::pair<size_t, uint32_t>> nbrs;
            nbrs.reserve(mol[atomIdx].degree());
            for (auto&& neighbor : mol[atomIdx].neighbors()) {
                auto bond = mol.bond(atomIdx, neighbor);

                roundAtomNeighborhoods[atomIdx][bond.index()] = 1;
                roundAtomNeighborhoods[atomIdx] |= atom_neighborhood[neighbor];

                size_t bt = 1;
                if (!useBondTypes) {
                    nbrs.push_back({0, invariants[neighbor]});
                    continue;
                }

                // REPLACE with switch statement
                if (!useChirality || bond.order() != Bond::DOUBLE /*|| bond.order() == Bond::STEREONONE*/) {
                    bt = static_cast<size_t>(bond.order());
                } else {
                    const size_t stereoOffset = 100;
                    const size_t bondTypeOffset = 10;
                    bt = stereoOffset +
                        bondTypeOffset * static_cast<size_t>(bond.order()) +
                        static_cast<size_t>(bond.order());
                }

                nbrs.push_back({bt, invariants[neighbor]});
            }

            // sort the neighbor list:
            std::sort(nbrs.begin(), nbrs.end());

            // and now calculate the new invariant and test if the atom is newly "chiral"
            auto invar = layer;
            boost::hash_combine(invar, invariants[atomIdx]);
            bool looksChiral = atom.chirality() != Atom::UNSPECIFIED;
            for (const auto& it : nbrs) {
                // add the contribution to the new invariant:
                boost::hash_combine(invar, it);

                // update our "chirality":
                if (useChirality && looksChiral && chiralAtoms[atomIdx]) {
                    if (it.first != static_cast<int32_t>(Bond::SINGLE)) {
                        looksChiral = false;
                    } /* else if (it != nbrs.begin() && it->second == (it - 1)->second) {
                        looksChiral = false;
                    } */
                }
            }

            if (useChirality && looksChiral) {
                chiralAtoms[atomIdx] = true;
                // add an extra value to the invariant to reflect chirality:
                if (atom.chirality() == Atom::CIP_R) {
                    boost::hash_combine(invar, 3);
                } else if (atom.chirality() == Atom::CIP_S) {
                    boost::hash_combine(invar, 2);
                } else {
                    boost::hash_combine(invar, 1);
                }
            }

            roundInvariants[atomIdx] = invar;
            round_neighborhoods.push_back(
                std::make_tuple(
                    roundAtomNeighborhoods[atomIdx],
                    roundInvariants[atomIdx],
                    atomIdx
                )
            );

            // we have seen this exact environment before, this atom
            // is now out of consideration:
            if (contains(neighborhoods, roundAtomNeighborhoods[atomIdx])) {
                deadAtoms[atomIdx] = 1;
            }
        }

        std::sort(round_neighborhoods.begin(), round_neighborhoods.end());
        for (auto& iter : round_neighborhoods) {
            // if we haven't seen this exact environment before, update the fingerprint:
            if (!contains(neighborhoods, std::get<0>(iter))) {
                if (onlyNonzeroInvariants && invariants_in[std::get<2>(iter)] == 0) {
                    continue;
                }
                if (include_atom[std::get<2>(iter)]) {
                    auto bit = res.increase(std::get<1>(iter));
                    atomsSettingBits[bit].push_back({std::get<2>(iter), layer + 1});
                }
                if (fromAtoms.empty() || contains(fromAtoms, std::get<2>(iter))) {
                    neighborhoods.push_back(std::get<0>(iter));
                }
            } else {
                // we have seen this exact environment before, this atom
                // is now out of consideration:
                deadAtoms[std::get<2>(iter)] = 1;
            }
        }

        // the invariants from this round become the global invariants:
        std::copy(roundInvariants.begin(), roundInvariants.end(),
                  invariants.begin());

        atom_neighborhood = roundAtomNeighborhoods;
    }

    return res;
}
}
