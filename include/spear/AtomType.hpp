// License:

#ifndef SPEAR_ATOMTYPE_HPP
#define SPEAR_ATOMTYPE_HPP

#include <vector>
#include <stdexcept>
#include <functional>   // std::bad_function_call

#include "spear/exports.hpp"
#include "spear/Constants.hpp"

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

class SPEAR_EXPORT AtomType : protected std::vector<size_t> {
protected:
    using super = std::vector<size_t>;
    void set(size_t idx, size_t val) {
        super::operator[](idx) = val;
    }
    const size_t& get(size_t idx) const {
        return super::operator[](idx);
    }
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

    /// Is the given atom aromatic in the typing scheme?
    virtual bool is_aromatic(size_t atom_id) const = 0;

    /// Is the given atom forced to be planar in the typing scheme?
    virtual bool is_planar(size_t atom_id) const = 0;

    /// Return the hybridization of the atom
    virtual Hybridization hybridization(size_t atom_id) const = 0;

    /// Add a new atom to the atom type vector
    virtual size_t add_atom(size_t new_idx) = 0;

    /// Retype the atoms idx1 and idx2 based on bo
    virtual void add_bond(size_t idx1, size_t idx2, Bond::Order bo) = 0;

    /// Remove an atom from the atom type vector
    virtual void remove_atom(size_t idx) = 0;

    /// Return a vector of all the atom type ids represented by the atom type
    const vector<size_t>& as_vec() const{
        return *this;
    }

    using super::size;

    /// Beginning of all atomtypes
    using super::cbegin;
    using super::cend;

    /// Retreive the type of a given typed atom
    using super::operator[];

    super::const_iterator begin() const {
        return cbegin();
    }

    super::const_iterator end() const {
        return cend();
    }

    friend bool operator==(const AtomType& at1, const AtomType& at2) {
        return at1.as_vec() == at2.as_vec();
    }
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
size_t atomtype_id_for_name(const std::string& name) {
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
