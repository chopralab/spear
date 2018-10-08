// License:

#ifndef SPEAR_ATOMTYPE_HPP
#define SPEAR_ATOMTYPE_HPP

#include "spear/exports.hpp"
#include "spear/Molecule.hpp"

namespace Spear {

class SPEAR_EXPORT AtomType {
public:
    /// Get the name of the set of atom types
    virtual std::string name() const = 0;

    /// Type the given molecule using 3D geometry
    virtual void type_atoms_3d(const Molecule& mol) = 0;

    /// Type the given molecule using bond orders
    virtual void type_atoms_order(const Molecule& mol) = 0;

    /// Get the name of a given type
    virtual std::string name(size_t id) const = 0;

    /// Get the id for a name
    virtual size_t id(const std::string& name) const = 0;

    /// Get all unique type ids
    virtual std::unordered_set<size_t> unique_ids() const = 0;

    /// Get all types
    virtual const std::vector<size_t>& all_ids() const = 0;
};

}

#endif
