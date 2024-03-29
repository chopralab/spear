// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_SYBYL_HPP
#define SPEAR_SYBYL_HPP

#include "spear/AtomType.hpp"

namespace Spear {

class AtomVertex;
class Molecule;

class SPEAR_EXPORT Sybyl : public AtomType {
public:

    Sybyl(const Molecule& mol, TypingMode mode);

    const std::string& name() const override {
        return name_;
    }

    bool is_aromatic(size_t atom_id) const override;

    bool is_planar(size_t atom_id) const override;

    size_t add_atom(size_t idx) override;

    void add_bond(size_t idx1, size_t idx2, Bond::Order bo) override;

    void remove_atom(size_t idx) override {
        erase(begin() + static_cast<std::ptrdiff_t>(idx));
    }

    Hybridization hybridization(size_t atom_id) const override;

private:

    void type_atoms_3d_();
    void type_atoms_topo_(const AtomVertex& atom);
    void type_atoms_topo_();

    size_t assign_carbon_topo_(const AtomVertex& atom);
    size_t assign_nitrogen_topo_(const AtomVertex& atom);

    size_t assign_carbon_3d_(const AtomVertex& atom);
    size_t assign_nitrogen_3d_(const AtomVertex& atom);
    size_t assign_oxygen_(const AtomVertex& atom);
    size_t assign_sulfur_(const AtomVertex& atom);

    /// Original molecule to be typed
    const Molecule& mol_;

    /// Generated name
    std::string name_;
};

template<> std::string SPEAR_EXPORT atomtype_name_for_id<Sybyl>(size_t id);

template<> size_t SPEAR_EXPORT atomtype_id_for_name<Sybyl>(const std::string& name);

template<> size_t SPEAR_EXPORT atomtype_id_count<Sybyl>();

}

#endif
