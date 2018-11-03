// License:

#ifndef SPEAR_ATOMTYPE_HPP
#define SPEAR_ATOMTYPE_HPP

#include <functional>   // std::bad_function_call

#include "spear/exports.hpp"
#include "spear/Molecule.hpp"

namespace Spear {

class SPEAR_EXPORT AtomType {
public:
    /// Type the given molecule using 3D geometry
    virtual std::vector<size_t> type_atoms_3d(const Molecule& mol) = 0;

    /// Type the given molecule using bond orders
    virtual std::vector<size_t> type_atoms_order(const Molecule& mol) = 0;
};

template<class Format>
std::string atomtype_name() {
    throw std::bad_function_call("atomtype_name is unimplemented for this format");
}

template<class Format>
std::string atomtype_name_for_id(size_t id) {
    throw std::bad_function_call("atomtype_name_for_id is unimplemented for this format");
}

template<class Format>
size_t atomtype_id_for_name(std::string name) {
    throw std::bad_function_call("atomtype_id_for_name is unimplemented for this format");
}

}

#endif
