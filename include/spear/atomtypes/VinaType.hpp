#ifndef SPEAR_VINATYPE_HPP
#define SPEAR_VINATYPE_HPP

#include "spear/AtomType.hpp"

namespace Spear {

class Molecule;

class SPEAR_EXPORT VinaType final : public AtomType {
public:

    enum TypingMode : size_t {
        ACCEPTOR_WATER = 1,
        DONAR_WATER = 2,
        HIS_NE2_PROTONATED = 4,
        HIS_ND1_PROTONATED = 8,
    };

    VinaType(const Molecule& mol, uint64_t mode = ACCEPTOR_WATER | DONAR_WATER | HIS_NE2_PROTONATED);

    const std::string& name() const override {
        return name_;
    }

    size_t add_atom(size_t new_idx) override;

    void add_bond(size_t idx1, size_t idx2, Bond::Order bo) override;

    void remove_atom(size_t idx) override {
        erase(this->begin() + static_cast<std::ptrdiff_t>(idx));
    }

    bool is_aromatic(size_t atom_id) const override {
        return false;
    }

    bool is_planar(size_t atom_id) const override {
        return false;
    }

    Hybridization hybridization(size_t atom_id) const override {
        return Hybridization::SP3;
    }

    static bool is_hydrophobic(size_t xs);
    static bool is_acceptor(size_t xs);
    static bool is_donor(size_t xs);
    static bool is_skip(size_t xs);

private:

    const Molecule& mol_;
    uint64_t mode_;

    const std::string name_ = "vina";
};

template<> std::string SPEAR_EXPORT atomtype_name_for_id<VinaType>(size_t id);

template<> size_t SPEAR_EXPORT atomtype_id_for_name<VinaType>(const std::string& name);

template<> size_t SPEAR_EXPORT atomtype_id_count<VinaType>();

template<> double SPEAR_EXPORT van_der_waals<VinaType>(size_t id);

}

#endif
