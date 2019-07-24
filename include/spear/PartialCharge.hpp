// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_PARTIALCHARGE_HPP
#define SPEAR_PARTIALCHARGE_HPP

#include <vector>
#include <stdexcept>
#include <functional>   // std::bad_function_call

#include "spear/exports.hpp"
#include "spear/Constants.hpp"

namespace Spear {

class SPEAR_EXPORT PartialCharge : protected std::vector<double> {
protected:
    using super = std::vector<double>;
    void set(size_t idx, double val) {
        super::operator[](idx) = val;
    }
    const double& get(size_t idx) const {
        return super::operator[](idx);
    }

    PartialCharge(super charges) : super(charges) {}
public:

    PartialCharge() {}

    virtual ~PartialCharge() = default;

    /// Get the name of partial charge. This maybe dependant on how the charges were
    /// initialized or how it is going to be used.
    virtual const std::string& name() const = 0;    

    /// Add a new atom to the atom type vector
    virtual double add_atom(size_t new_idx) = 0;

    /// Remove an atom from the atom type vector
    virtual void remove_atom(size_t idx) = 0;

    /// Return a vector of all the atom type ids represented by the atom type
    const std::vector<double>& as_vec() const{
        return *this;
    }

    using super::size;

    /// Beginning of all atomtypes
    using super::cbegin;
    using super::cend;

    /// Retreive the type of a given charged atom
    using super::operator[];

    super::const_iterator begin() const {
        return cbegin();
    }

    super::const_iterator end() const {
        return cend();
    }

    friend bool operator==(const PartialCharge& at1, const PartialCharge& at2) {
        return at1.as_vec() == at2.as_vec();
    }
};

}

#endif
