// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_GAFF_HPP
#define SPEAR_GAFF_HPP

#include "spear/AtomType.hpp"

#include <map>

namespace Spear {

class Molecule;

class SPEAR_EXPORT GAFF : public AtomType {
public:

    GAFF(const Molecule& mol);

    const std::string& name() const override {
        return name_;
    }

    size_t add_atom(size_t new_idx) override;

    void add_bond(size_t idx1, size_t idx2, Bond::Order bo) override;

    void remove_atom(size_t idx) override {
        erase(begin() + static_cast<std::ptrdiff_t>(idx));
    }

    bool is_aromatic(size_t atom_id) const override;

    bool is_planar(size_t atom_id) const override;

    Hybridization hybridization(size_t /*atom_id*/) const override {
		return Hybridization::FORCED;
    }

private:

    size_t add_atom_(size_t new_idx);

    const Molecule& mol_;

    const std::string name_ = "gaff";
};

template<> std::string SPEAR_EXPORT atomtype_name_for_id<GAFF>(size_t id);

template<> size_t SPEAR_EXPORT atomtype_id_for_name<GAFF>(const std::string& name);

template<> size_t SPEAR_EXPORT atomtype_id_count<GAFF>();

template<> double SPEAR_EXPORT van_der_waals<GAFF>(size_t id);

}

#endif
