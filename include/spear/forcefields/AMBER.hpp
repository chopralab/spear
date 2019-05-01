// License:

#ifndef SPEAR_AMBER_HPP
#define SPEAR_AMBER_HPP

#include "spear/Forcefield.hpp"

#include <unordered_map>

namespace pugi {
    class xml_node;
}

namespace chemfiles {
    class Residue;
}

namespace Spear {

class AtomVertex;

class SPEAR_EXPORT AMBER : public BondedForcefield, public NonBondedForcefield {
public:

    struct AtomType {
        size_t atom_class, atom_template;
        double mass, charge, sigma, epsilon;
    };

    using Residue = std::unordered_map<std::string, AtomType>;

    struct Bond {
        double length, k;
    };

    struct Angle {
        double theta, k;
    };

    struct PeriodicForce {
        double phase, k, periodicity;
    };

    using Torsion = std::vector<PeriodicForce>;

    AMBER(std::istream& input);

    virtual ~AMBER(){}

    void add_xml_file(std::istream& input) {
        read_dat_file_(input);
    }

    void add_forces(const std::vector<const Molecule*>& mols, OpenMM::System& system) const override;

    void add_forces(const Molecule& mol, OpenMM::System& system) const override;

    std::vector<double> masses(const Molecule& mol) const override;

    void set_cutoff(double cutoff);

    void set_method(NonbondedMethod method) {
        method_ = method;
    }

    void reset() {
        non_bond_force_ = -1;
        bond_force_ = -1;
        angle_force_ = -1;
        torsion_force_ = -1;
        improper_force_ = -1;
    }

private:

    template<size_t N, typename T>
    bool get_residue_pair(const Molecule& mol,
                          const T& atoms, bool use_class,
                          std::array<size_t, N>& out) const;

    template<typename Func>
    void apply_function_to_atoms_in_residue(const Molecule& mol, Func&& func) const;

    const AtomType& find_atom_type(const AtomVertex& atom, const chemfiles::Residue& res) const;

    size_t get_class_(const std::string& s) const;
    size_t get_type_(const std::string& s) const;

    template<size_t N>
    bool read_classes_(pugi::xml_node& node, const std::string& prefix,
                       std::array<size_t, N>& ids, size_t offset = 0) const;

    std::unordered_map<std::string, size_t> class_map_;
    std::unordered_map<std::string, size_t> type_map_;
    std::vector<std::string> class_names_;
    std::vector<std::string> type_names_;
    std::vector<AtomType> type_templates_;

    void read_dat_file_(std::istream& input);

    void read_atom_types_(pugi::xml_node& input);
    std::unordered_map<size_t, size_t> class_to_type_;
    std::unordered_map<size_t, size_t> type_to_class_;

    void read_residues_(pugi::xml_node& input);
    using ResidueMap = std::unordered_map<std::string, Residue>;
    ResidueMap residues_;

    void read_bonds_(pugi::xml_node& input);
    std::unordered_map<bond_type, Bond, bond_type_hash, bond_type_equal> type_to_bond_;

    void read_angles_(pugi::xml_node& input);
    std::unordered_map<angle_type, Angle, angle_type_hash, angle_type_equal> type_to_angle_;

    void read_torsions_(pugi::xml_node& input);
    std::unordered_map<torsion_type, Torsion, torsion_type_hash, torsion_type_equal> type_to_torsion_;
    std::unordered_map<improper_type, Torsion, improper_type_hash, improper_type_equal> type_to_improper_;

    mutable int non_bond_force_ = -1;
    mutable int bond_force_ = -1;
    mutable int angle_force_ = -1;
    mutable int torsion_force_ = -1;
    mutable int improper_force_ = -1;

    size_t class_ids = 1000001;
    size_t type_ids = 1;

    double coulombic14scale_;
    double lj14scale_;

    bool res_charge_ = false;
    bool res_sigma_ = false;
    bool res_epsilon_ = false;

    double cutoff_;
    NonbondedMethod method_;
};

}

#endif
