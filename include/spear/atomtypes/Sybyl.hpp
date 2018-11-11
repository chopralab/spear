#ifndef SPEAR_SYBYL_HPP
#define SPEAR_SYBYL_HPP

#include "spear/AtomType.hpp"

namespace Spear {

class SPEAR_EXPORT Sybyl : public AtomType {
public:

    Sybyl(const Molecule& mol);

    void type_atoms_3d() override;
    void type_atoms_order() override;

    const std::vector<size_t>& all_types() const override {
        return atom_types_;
    }

    size_t operator[](size_t atom_id) const override {
        return atom_types_[atom_id];
    }

    std::vector<size_t>::const_iterator cbegin() const override {
        return atom_types_.cbegin();
    }

    std::vector<size_t>::const_iterator cend() const override {
        return atom_types_.cend();
    }

private:
    size_t assign_carbon_ord_(AtomVertex atom);
    size_t assign_nitrogen_ord_(AtomVertex atom);

    size_t assign_carbon_3d_(AtomVertex atom);
    size_t assign_nitrogen_3d_(AtomVertex atom);
    size_t assign_oxygen_(AtomVertex atom);
    size_t assign_sulfur_(AtomVertex atom);

    /// Vector to hold all atom types
    std::vector<size_t> atom_types_;

    /// Original molecule to be typed
    const Molecule& mol_;
};

template<> std::string atomtype_name<Sybyl>();

template<> std::string atomtype_name_for_id<Sybyl>(size_t id);

template<> size_t atomtype_id_for_name<Sybyl>(std::string name);

template<> size_t atomtype_id_count<Sybyl>();

}

#endif
