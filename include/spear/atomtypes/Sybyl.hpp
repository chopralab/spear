#ifndef SPEAR_SYBYL_HPP
#define SPEAR_SYBYL_HPP

#include "spear/AtomType.hpp"

namespace Spear {

class SPEAR_EXPORT Sybyl : public AtomType {
public:

    std::string name() const override;

    void type_atoms_3d(const Molecule& mol) override;

    void type_atoms_order(const Molecule& mol) override;

    std::string name(size_t id) const override;

    size_t id(const std::string& name) const override;

    std::unordered_set<size_t> unique_ids() const override;

    const std::vector<size_t>& all_ids() const override;

private:
    std::vector<size_t> atom_types_;

    size_t assign_carbon_ord_(const Molecule& mol, AtomVertex atom);
    size_t assign_nitrogen_ord_(const Molecule& mol, AtomVertex atom);

    size_t assign_carbon_3d_(const Molecule& mol, AtomVertex atom);
    size_t assign_nitrogen_3d_(const Molecule& mol, AtomVertex atom);
    size_t assign_oxygen_(const Molecule& mol, AtomVertex atom);
    size_t assign_sulfur_(const Molecule& mol, AtomVertex atom);
};

}

#endif
