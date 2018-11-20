// License:

#ifndef SPEAR_ATOMTYPE_HPP
#define SPEAR_ATOMTYPE_HPP

#include <vector>
#include <functional>   // std::bad_function_call

#include "spear/exports.hpp"

namespace Spear {

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

    /// Retreive the type of a given typed atom
    virtual size_t operator[](size_t atom_id) const = 0;

    /// Beginning of all atomtypes
    virtual std::vector<size_t>::const_iterator cbegin() const = 0;

    /// Ending of all atomtypes
    virtual std::vector<size_t>::const_iterator cend() const = 0;
};

template<class Format>
std::string atomtype_name_for_id(size_t id) {
    throw std::bad_function_call("atomtype_name_for_id is unimplemented for this format");
}

template<class Format>
size_t atomtype_id_for_name(std::string name) {
    throw std::bad_function_call("atomtype_id_for_name is unimplemented for this format");
}

template<class Format>
size_t atomtype_id_count() {
    throw std::bad_function_call("atomtype_id_count is unimplemented for this format");
}

template<class Format>
double van_der_waals(size_t id) {
    throw std::bad_function_call("van_der_waals is unimplemented for this format");
}

}

#endif
