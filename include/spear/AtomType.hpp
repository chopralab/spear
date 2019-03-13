// License:

#ifndef SPEAR_ATOMTYPE_HPP
#define SPEAR_ATOMTYPE_HPP

#include <vector>
#include <functional>   // std::bad_function_call

#include "spear/exports.hpp"

namespace Spear {

enum Hybridization : uint64_t {
    UNKNOWN = 0,
    FORCED  = 1,
    SP      = 2,
    SP2     = 3,
    SP3     = 4,
    SP3D    = 5,
    SP3D2   = 6,
};

class SPEAR_EXPORT AtomType {
public:

    enum TypingMode {
        GEOMETRY,
        TOPOLOGY,
        AUTO
    };

    virtual ~AtomType() = default;

    /// Get the name of atom type. This maybe dependant on how the type was
    /// initialized or how it is going to be used.
    virtual const std::string& name() const = 0;    

    /// Return a vector of all the atom type ids represented by the atom type
    virtual const std::vector<size_t>& all_types() const = 0;

    /// Is the given atom aromatic in the typing scheme?
    virtual bool is_aromatic(size_t atom_id) const = 0;

    /// Return the hybridization of the atom
    virtual Hybridization hybridization(size_t atom_id) const = 0;

    /// Add a new atom to the atom type vector
    virtual size_t add_atom(size_t new_idx) = 0;

    /// Remove an atom from the atom type vector
    virtual void remove_atom(size_t idx) = 0;

    /// Retreive the type of a given typed atom
    virtual size_t operator[](size_t atom_id) const = 0;

    /// Beginning of all atomtypes
    virtual std::vector<size_t>::const_iterator cbegin() const = 0;

    /// Ending of all atomtypes
    virtual std::vector<size_t>::const_iterator cend() const = 0;
};

class FormatFeatureUnimplemented : public std::logic_error {
public:
    FormatFeatureUnimplemented(const std::string& s) :
        std::logic_error(s + " is unimplemented for this format.") {}
};

template<class Format>
std::string atomtype_name_for_id(size_t id) {
    throw FormatFeatureUnimplemented("atomtype_name_for_id");
}

template<class Format>
size_t atomtype_id_for_name(std::string name) {
    throw FormatFeatureUnimplemented("atomtype_id_for_name");
}

template<class Format>
size_t atomtype_id_count() {
    throw FormatFeatureUnimplemented("atomtype_id_count");
}

template<class Format>
double van_der_waals(size_t id) {
    throw FormatFeatureUnimplemented("van_der_waals");
}

}

#endif
