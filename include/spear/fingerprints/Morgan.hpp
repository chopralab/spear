#include <cctype>
#include <vector>

#include "spear/Fingerprint.hpp"
#include "spear/exports.hpp"

namespace Spear{

class Molecule;

Fingerprint SPEAR_EXPORT calcFingerprint(
    const Molecule& mol, std::uint64_t radius,
    const std::vector<size_t>& invariants_in = std::vector<size_t>(),
    const std::vector<size_t>& fromAtoms = std::vector<size_t>(),
    size_t finger_print_length = 2048,
    bool useChirality = false, bool useBondTypes = true,
    bool onlyNonzeroInvariants = false
);

}
