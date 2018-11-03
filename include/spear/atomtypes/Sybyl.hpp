#ifndef SPEAR_SYBYL_HPP
#define SPEAR_SYBYL_HPP

#include "spear/AtomType.hpp"

namespace Spear {

class SPEAR_EXPORT Sybyl : public AtomType {
public:

    std::vector<size_t> type_atoms_3d(const Molecule& mol) override;

    std::vector<size_t> type_atoms_order(const Molecule& mol) override;

private:
    size_t assign_carbon_ord_(const Molecule& mol, AtomVertex atom);
    size_t assign_nitrogen_ord_(const Molecule& mol, AtomVertex atom);

    size_t assign_carbon_3d_(const Molecule& mol, AtomVertex atom);
    size_t assign_nitrogen_3d_(const Molecule& mol, AtomVertex atom);
    size_t assign_oxygen_(const Molecule& mol, AtomVertex atom);
    size_t assign_sulfur_(const Molecule& mol, AtomVertex atom);
};

template<> std::string atomtype_name<Sybyl>();

template<> std::string atomtype_name_for_id<Sybyl>(size_t id);

template<> size_t atomtype_id_for_name<Sybyl>(std::string name);

}

#endif
