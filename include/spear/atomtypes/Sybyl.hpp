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

    const std::vector<size_t>& all_types() const override {
        return atom_types_;
    }

    bool is_aromatic(size_t atom_id) const override;

    bool is_planar(size_t atom_id) const override;

    size_t add_atom(size_t idx) override;

    void remove_atom(size_t idx) {
        atom_types_.erase(atom_types_.begin() + idx);
    }

    Hybridization hybridization(size_t atom_id) const override;

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

    void type_atoms_3d_();
    void type_atoms_topo_();

    size_t assign_carbon_topo_(AtomVertex& atom);
    size_t assign_nitrogen_topo_(AtomVertex& atom);

    size_t assign_carbon_3d_(AtomVertex& atom);
    size_t assign_nitrogen_3d_(AtomVertex& atom);
    size_t assign_oxygen_(AtomVertex& atom);
    size_t assign_sulfur_(AtomVertex& atom);

    /// Vector to hold all atom types
    std::vector<size_t> atom_types_;

    /// Original molecule to be typed
    const Molecule& mol_;

    /// Generated name
    std::string name_;
};

template<> std::string SPEAR_EXPORT atomtype_name_for_id<Sybyl>(size_t id);

template<> size_t SPEAR_EXPORT atomtype_id_for_name<Sybyl>(std::string name);

template<> size_t SPEAR_EXPORT atomtype_id_count<Sybyl>();

}

#endif
