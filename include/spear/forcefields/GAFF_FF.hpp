// License:

#ifndef SPEAR_GAFF_FF_HPP
#define SPEAR_GAFF_FF_HPP

#include "spear/Forcefield.hpp"

#include <map>

namespace Spear {

class SPEAR_EXPORT GAFF_FF : public BondedForcefield {
public:

    struct Atom {
        double mass, polarizability, sigma, epsilon;
    };

    struct Bond {
        double length, k;
    };

    struct Angle {
        double theta, k;
    };

    struct Torsion {
        double phase, k, periodicity;
    };

    struct Improper {
        double phase, k, periodicity;
    };

    GAFF_FF(std::istream& input);

    virtual ~GAFF_FF(){}

    void add_forces(const Molecule& mol, OpenMM::System& system) const override;

    std::vector<double> masses(const Molecule& mol) const override;

    size_t num_atom_types() const {
        return type_to_atom_.size();
    }

    size_t num_bonds() const {
        return type_to_bond_.size();
    }

    size_t num_angles() const {
        return type_to_angle_.size();
    }

    size_t num_torsions() const {
        return type_to_torsion_.size();
    }

    size_t num_impropers() const {
        return type_to_improper_.size();
    }

    void reset() {
        non_bond_force_ = -1;
        bond_force_ = -1;
        angle_force_ = -1;
        torsion_force_ = -1;
        improper_force_ = -1;
    }

private:
    void read_dat_file_(std::istream& input);

    void read_atom_types_(std::istream& input);
    std::unordered_map<size_t, Atom> type_to_atom_;

    void read_bonds_(std::istream& input);
    std::unordered_map<bond_type, Bond, bond_type_hash, bond_type_equal> type_to_bond_;

    void read_angles_(std::istream& input);
    std::unordered_map<angle_type, Angle, angle_type_hash, angle_type_equal> type_to_angle_;

    void read_torsions_(std::istream& input);
    std::unordered_map<torsion_type, Torsion, torsion_type_hash, torsion_type_equal> type_to_torsion_;

    void read_impropers_(std::istream& input);
    std::unordered_map<improper_type, Torsion, improper_type_hash, improper_type_equal> type_to_improper_;

    void read_lj_(std::istream& input);

    mutable int non_bond_force_ = -1;
    mutable int bond_force_ = -1;
    mutable int angle_force_ = -1;
    mutable int torsion_force_ = -1;
    mutable int improper_force_ = -1;
};

}

#endif
