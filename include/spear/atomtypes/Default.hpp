// License:

#ifndef SPEAR_DEFAULT_HPP
#define SPEAR_DEFAULT_HPP

#include "spear/AtomType.hpp"

#include <map>

namespace Spear {

class Molecule;

class SPEAR_EXPORT Default final : public AtomType {
public:

    Default(const Molecule& mol);

    const std::string& name() const override {
        return name_;
    }

    size_t add_atom(size_t new_idx) override;

    void add_bond(size_t idx1, size_t idx2, Bond::Order bo) override;

    void remove_atom(size_t idx) override {
        erase(this->begin() + static_cast<std::ptrdiff_t>(idx));
        hybridizations_.erase(hybridizations_.begin() + static_cast<std::ptrdiff_t>(idx));
    }

    bool is_aromatic(size_t atom_id) const override;

    bool is_planar(size_t atom_id) const override;

    Hybridization hybridization(size_t atom_id) const override {
        return hybridizations_[atom_id];
    }

private:

    const Molecule& mol_;

    /// number of heavy atoms bonded
    std::vector<Hybridization> hybridizations_;

    const std::string name_ = "default";
};

template<> std::string SPEAR_EXPORT atomtype_name_for_id<Default>(size_t id);

template<> size_t SPEAR_EXPORT atomtype_id_for_name<Default>(const std::string& name);

template<> size_t SPEAR_EXPORT atomtype_id_count<Default>();

template<> double SPEAR_EXPORT van_der_waals<Default>(size_t id);

}

#endif
