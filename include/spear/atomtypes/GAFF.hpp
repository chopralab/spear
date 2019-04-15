// License:

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

    const std::vector<size_t>& all_types() const override {
        return atom_types_;
    }

    size_t add_atom(size_t new_idx) override;

    void add_bond(size_t idx1, size_t idx2, Bond::Order bo) override;

    void remove_atom(size_t idx) override {
        atom_types_.erase(atom_types_.begin() + idx);
    }

    size_t operator[](size_t atom_id) const override {
        return atom_types_[atom_id];
    }

    bool is_aromatic(size_t atom_id) const override;

    bool is_planar(size_t atom_id) const override;

    Hybridization hybridization(size_t atom_id) const override {
		return Hybridization::FORCED;
    }

    std::vector<size_t>::const_iterator cbegin() const override {
        return atom_types_.cbegin();
    }

    std::vector<size_t>::const_iterator cend() const override {
        return atom_types_.cend();
    }

private:

    const Molecule& mol_;
    std::vector<size_t> atom_types_;

    const std::string name_ = "gaff";
};

template<> std::string SPEAR_EXPORT atomtype_name_for_id<GAFF>(size_t id);

template<> size_t SPEAR_EXPORT atomtype_id_for_name<GAFF>(std::string name);

template<> size_t SPEAR_EXPORT atomtype_id_count<GAFF>();

template<> double SPEAR_EXPORT van_der_waals<GAFF>(size_t id);

}

#endif
