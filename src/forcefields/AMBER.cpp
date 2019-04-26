#include "spear/forcefields/AMBER.hpp"
#include "spear/Molecule.hpp"
#include "spear/Molecule_impl.hpp"
#include "spear/Geometry.hpp"

#include "OpenMM.h"

#include "pugixml.hpp"

#include <string>
#include <iostream>

using namespace Spear;

template<typename T>
T& get_force(OpenMM::System& system, int& force_id) {
    T* force = nullptr;
    if (force_id == -1) {
        force = new T();
        force_id = system.addForce(force);
    } else {
        force = dynamic_cast<T*>(&system.getForce(force_id));
    }

    return *force;
}

template<typename Func>
void AMBER::apply_function_to_atoms_in_residue(const Molecule& mol, Func&& func) const {
    auto& topo = mol.topology();
    for (auto& res : topo.residues()) {
        bool is_nterminal = res.get<chemfiles::Property::BOOL>("is_n_terminal").value_or(false);
        bool is_cterminal = res.get<chemfiles::Property::BOOL>("is_c_terminal").value_or(false);
        auto name = res.name();
        if (is_nterminal) {
            name = "N" + name;
        } else if (is_cterminal) {
            name = "C" + name;
        }

        auto res_iter = residues_.find(name);
        if (res_iter == residues_.end()) {
            std::cerr << "Resiude not found: " << name << std::endl;
            continue;
        }

        for (auto& atom : res) {
            const auto& atom_class = mol[atom];
            auto atom_iter = res_iter->second.find(atom_class.name());
            if (atom_iter == res_iter->second.end()) {
                std::cerr << "Atom not found: " << atom_class.name()
                          << " in " << res.name() << std::endl;
                continue;
            }
            func(atom, atom_iter->second);
        }
    }
}

template<size_t N, typename T>
bool AMBER::get_residue_pair(const Molecule& mol,
                             const T& atoms, bool use_class,
                             std::array<size_t, N>& out) const {

    for (size_t i = 0; i < N; ++i) {
        auto res = mol.topology().residue_for_atom(atoms[i]);
        if (!res) {
            return false;
        }
        out[i] = use_class ? find_atom_type(mol[atoms[i]], *res).atom_class
                           : find_atom_type(mol[atoms[i]], *res).atom_template;
    }

    return true;
}

const AMBER::AtomType& AMBER::find_atom_type(const AtomVertex& atom,
                                             const chemfiles::Residue& res) const {

    bool is_nterminal = res.get<chemfiles::Property::BOOL>("is_n_terminal").value_or(false);
    bool is_cterminal = res.get<chemfiles::Property::BOOL>("is_c_terminal").value_or(false);
    auto name = res.name();
    if (is_nterminal) {
        name = "N" + name;
    } else if (is_cterminal) {
        name = "C" + name;
    }
    auto res_iter = residues_.find(name);
    if (res_iter == residues_.end()) {
        throw std::runtime_error("Resiude not found: " + name);
    }

    auto atom_iter = res_iter->second.find(atom.name());
    if (atom_iter == res_iter->second.end()) {
        throw std::runtime_error("Atom not found: " + atom.name());
    }

    return atom_iter->second;
}


AMBER::AMBER(std::istream& input) {
    read_dat_file_(input);
}

void AMBER::add_forces(const Molecule& mol, OpenMM::System& system) const {
    if (static_cast<int>(mol.size()) > system.getNumParticles()) {
        throw std::runtime_error("You must add all particles to the system first.");
    }

    auto added = static_cast<size_t>(system.getNumParticles()) - mol.size();

    auto& nonbond = get_force<OpenMM::NonbondedForce>(system, non_bond_force_);

    for (auto junk : mol) {
        nonbond.addParticle(0, 0, 0);
    }

    apply_function_to_atoms_in_residue(mol,
        [this, &nonbond](size_t atom, const AtomType& type){
            nonbond.setParticleParameters(static_cast<int>(atom), type.charge, type.sigma, type.epsilon);
        }
    );

    auto& topo = mol.topology();

    auto& hbond = get_force<OpenMM::HarmonicBondForce>(system, bond_force_);
    std::vector<std::pair<int,int>> bonds;
    for (auto bond : topo.bonds()) {
        bonds.emplace_back(static_cast<int>(bond[0] + added),
                           static_cast<int>(bond[1] + added)
        );

        std::array<size_t,2> class_;
        if (!get_residue_pair(mol, bond, false, class_)) {
            continue;
        }

        auto lookup = type_to_bond_.find(class_);
        if (lookup == type_to_bond_.end()) {
            get_residue_pair(mol, bond, true, class_);
            lookup = type_to_bond_.find(class_);
        }

        if (lookup == type_to_bond_.end()) {
            std::cerr << "Warning: bond '"
                      << mol[bond[0]].name() << " to "
                      << mol[bond[1]].name()
                      << " not found.\n";
            continue;
        }

        hbond.addBond(static_cast<int>(bond[0] + added),
                      static_cast<int>(bond[1] + added),
            lookup->second.length, lookup->second.k
        );
    }

    nonbond.createExceptionsFromBonds(bonds, coulombic14scale_, lj14scale_);

    auto& hangle = get_force<OpenMM::HarmonicAngleForce>(system, angle_force_);
    for (auto angle : mol.topology().angles()) {
        std::array<size_t, 3> class_;
        if (!get_residue_pair(mol, angle, false, class_)) {
            continue;
        }

        auto lookup = type_to_angle_.find(class_);
        if (lookup == type_to_angle_.end()) {
            get_residue_pair(mol, angle, true, class_);
            lookup = type_to_angle_.find(class_);
        }

        if (lookup == type_to_angle_.end()) {
            std::cerr << "Warning: angle '"
                      << " not found.\n";
            continue;
        }

        hangle.addAngle(static_cast<int>(angle[0] + added),
                        static_cast<int>(angle[1] + added),
                        static_cast<int>(angle[2] + added),
            lookup->second.theta, lookup->second.k
        );
    }

    auto& torsion = get_force<OpenMM::PeriodicTorsionForce>(system, torsion_force_);
    for (auto dihedral : mol.topology().dihedrals()) {
        std::array<size_t,4> class_;
        if (!get_residue_pair(mol, dihedral, false, class_)) {
            continue;
        }

        auto lookup = type_to_torsion_.find(class_);
        if (lookup == type_to_torsion_.end()) {
            lookup = type_to_torsion_.find({0, class_[1], class_[2], 0});
        }

        if (lookup == type_to_torsion_.end()) {
            get_residue_pair(mol, dihedral, true, class_);
            lookup = type_to_torsion_.find(class_);
        }
        if (lookup == type_to_torsion_.end()) {
            lookup = type_to_torsion_.find({0, class_[1], class_[2], 0});
        }

        if (lookup == type_to_torsion_.end()) {
            std::cerr << "Warning: torison '"
                      << " not found.\n";
            continue;
        }

        for (auto force : lookup->second) {

            torsion.addTorsion(static_cast<int>(dihedral[0] + added),
                               static_cast<int>(dihedral[1] + added),
                               static_cast<int>(dihedral[2] + added),
                               static_cast<int>(dihedral[3] + added),
                               static_cast<int>(std::abs(force.periodicity)),
                               force.phase,
                               force.k
            );
        }
    }

    for (auto improper : mol.topology().impropers()) {

        std::array<size_t,4> class_;
        if (!get_residue_pair(mol, improper, false, class_)) {
            continue;
        }

        auto lookup = type_to_improper_.find(class_);
        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[0], class_[2]});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[0], class_[3]});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[2], class_[3]});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[0], 0});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[2], 0});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[3], 0});
        }


        if (lookup == type_to_improper_.end()) {
            get_residue_pair(mol, improper, true, class_);
            lookup = type_to_torsion_.find(class_);
        }
        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[0], 0});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[2], 0});
        }

        if (lookup == type_to_improper_.end()) {
            lookup = type_to_improper_.find({0, class_[1], class_[3], 0});
        }
        // No warning, these are rare
        if (lookup == type_to_improper_.end()) {
            continue;
        }

        // Note the switch!
        for (auto force : lookup->second) {
            torsion.addTorsion(static_cast<int>(improper[1] + added), // <------- OpenMM swap
                               static_cast<int>(improper[0] + added),
                               static_cast<int>(improper[2] + added),
                               static_cast<int>(improper[3] + added),
                               static_cast<int>(std::abs(force.periodicity)),
                               force.phase,
                               force.k
            );
        }
    }
}

std::vector<double> AMBER::masses(const Molecule& mol) const {
    std::vector<double> result(mol.size());

    apply_function_to_atoms_in_residue(mol,
        [this, &result](size_t atom, const AtomType& type){
        result[atom] = type.mass;
    });

    return result;
}

size_t AMBER::get_class_(const std::string& s) const {
    if (s == "") {
        return 0; // wild card
    }

    auto iter = class_map_.find(s);
    if (iter == class_map_.end()) {
        throw std::runtime_error("Atom class " + s + " not found!");
    }

    return iter->second;
}

size_t AMBER::get_type_(const std::string& s) const {
    if (s == "") {
        return 0; // wild card
    }

    auto iter = type_map_.find(s);
    if (iter == type_map_.end()) {
        throw std::runtime_error("Atom class " + s + " not found!");
    }

    return iter->second;
}

template<size_t N>
bool AMBER::read_classes_(pugi::xml_node& node, const std::string& prefix,
                   std::array<size_t, N>& ids, size_t offset) const {
    for (size_t i = 0; i < N; ++i) {
        auto name = prefix + std::to_string(i + 1);
        auto attr = node.attribute(name.c_str());
        if (!attr) {
            return false;
        }
        if (prefix == "type") {
            ids[i] = get_type_(attr.as_string()) + offset;
        } else {
            ids[i] = get_class_(attr.as_string()) + offset;
        }
    }

    return true;
}

void AMBER::read_dat_file_(std::istream& input) {
    pugi::xml_document doc;
    auto result = doc.load(input);

    if (!result) {
        throw std::runtime_error("Error parsing AMBER XML");
    }

    auto ff = doc.child("ForceField");
    if (!ff) {
        throw std::runtime_error("XML File must have ForceField child");
    }

    read_atom_types_(ff);
    read_residues_(ff);
    read_bonds_(ff);
    read_angles_(ff);
    read_torsions_(ff);
}

void AMBER::read_atom_types_(pugi::xml_node& input) {
    auto at_node = input.child("AtomTypes");
    if (!at_node) {
        return;
    }

    for (auto type : at_node.children("Type")) {
        auto name = type.attribute("name").as_string();
        auto class_ = type.attribute("class").as_string();
        auto mass = type.attribute("mass").as_double();

        type_names_.push_back(name);
        type_map_.insert({name, type_ids - 1});

        auto class_iter = class_map_.find(class_);
        if (class_iter == class_map_.end()) {
            class_iter = class_map_.insert({class_, class_ids}).first;
            class_ids++;
            class_names_.push_back(class_);
        }

        class_to_type_.insert({class_iter->second, type_ids});
        type_to_class_.insert({type_ids, class_iter->second});
        type_templates_.push_back({class_iter->second, type_ids, mass, 0, 0, 0});
        type_ids++;
    }

    at_node = input.child("NonbondedForce");

    coulombic14scale_ = at_node.attribute("coulombic14scale").as_double(0.8333333333333334);
    lj14scale_ = at_node.attribute("coulombic14scale").as_double(0.5);

    for (auto use_attribute : at_node.children("UseAttributeFromResidue")) {
        std::string name = use_attribute.attribute("name").as_string();
        if (name == "charge") {
            res_charge_ = true;
        }
        if (name == "sigma") {
            res_sigma_ = true;
        }
        if (name == "epsilon") {
            res_epsilon_ = true;
        }
    }

    for (auto atom : at_node.children("Atom")) {
        auto type = get_type_(atom.attribute("type").as_string());
        
        if (!res_charge_) {
            type_templates_[type].charge = atom.attribute("charge").as_double();
        }

        if (!res_sigma_) {
            type_templates_[type].sigma = atom.attribute("sigma").as_double();
        }

        if (!res_epsilon_) {
            type_templates_[type].epsilon = atom.attribute("epsilon").as_double();
        }
    }
}

void AMBER::read_residues_(pugi::xml_node& input) {
    auto at_node = input.child("Residues");
    if (!at_node) {
        return;
    }

    for (auto res : at_node.children("Residue")) {
        auto rname = res.attribute("name").as_string();

        Residue new_res;
        for (auto atom : res.children("Atom")) {
            auto name = atom.attribute("name").as_string();
            auto type = get_type_(atom.attribute("type").as_string());

            auto new_type = new_res.insert({name, type_templates_[type]});

            new_type.first->second.atom_template = type;

            if (res_charge_) {
                new_type.first->second.charge = atom.attribute("charge").as_double();
            }
            if (res_sigma_) {
                new_type.first->second.sigma = atom.attribute("sigma").as_double();
            }
            if (res_epsilon_) {
                new_type.first->second.epsilon = atom.attribute("epsilon").as_double();
            }
        }
        residues_.emplace(rname, new_res);
    }
}

void AMBER::read_bonds_(pugi::xml_node& input) {
    auto at_node = input.child("HarmonicBondForce");
    if (!at_node) {
        return;
    }

    for (auto bond : at_node.children("Bond")) {
        try {
            auto length = bond.attribute("length").as_double();
            auto k = bond.attribute("k").as_double();

            std::array<size_t, 2> bond_class;
            if (read_classes_(bond, "type", bond_class)) {
                type_to_bond_.insert({bond_class, {length, k}});
            }

            if (read_classes_(bond, "class", bond_class)) {
                type_to_bond_.insert({bond_class, {length, k}});
            }            
        } catch (const std::exception& e) {
            std::cerr << "Warning: " << e.what() << std::endl;
        } catch (...) {
            std::cout << "unknown error" << std::endl;
        }
    }
}

void AMBER::read_angles_(pugi::xml_node& input) {
    auto at_node = input.child("HarmonicAngleForce");
    if (!at_node) {
        return;
    }

    for (auto angle : at_node.children("Angle")) {
        try {
            auto theta = angle.attribute("angle").as_double();
            auto k = angle.attribute("k").as_double();

            std::array<size_t, 3> angle_class;
            if (read_classes_(angle, "type", angle_class)) {
                type_to_angle_.insert({angle_class, {theta, k}});
            }

            if (read_classes_(angle, "class", angle_class)) {
                type_to_angle_.insert({angle_class, {theta, k}});
            }
        } catch(const std::exception& e) {
            std::cerr << "Warning: " << e.what() << std::endl;
        }
    }
}

std::vector<AMBER::PeriodicForce> read_periodic_force(const pugi::xml_node& input) {
    size_t id = 1;
    std::vector<AMBER::PeriodicForce> forces;
    do {
        std::string periodicity = "periodicity";
        std::string phase = "phase";
        std::string k = "k";

        periodicity += std::to_string(id);
        phase += std::to_string(id);
        k += std::to_string(id);

        auto p_node = input.attribute(periodicity.c_str());
        auto f_node = input.attribute(phase.c_str());
        auto k_node = input.attribute(k.c_str());

        if (!p_node || !f_node || !k_node) {
            break;
        }

        forces.push_back({f_node.as_double(),
                          k_node.as_double(),
                          p_node.as_double()}
        );

        id++;
    } while(true);

    return forces;
}

void AMBER::read_torsions_(pugi::xml_node& input) {
    auto at_node = input.child("PeriodicTorsionForce");
    if (!at_node) {
        return;
    }

    for (auto proper : at_node.children("Proper")) {
        try {
            auto forces = read_periodic_force(proper);

            std::array<size_t, 4> torsion_class;
            if (read_classes_(proper, "type", torsion_class)) {
                type_to_torsion_.insert({torsion_class, std::move(forces)});
                continue;
            }

            if (read_classes_(proper, "class", torsion_class)) {
                type_to_torsion_.insert({torsion_class, std::move(forces)});
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: " << e.what() << std::endl;
        }
    }

    for (auto improper : at_node.children("Improper")) {
        try {
            auto forces = read_periodic_force(improper);

            std::array<size_t, 4> torsion_class;
            if (read_classes_(improper, "type", torsion_class)) {
                std::swap(torsion_class[0], torsion_class[1]);
                type_to_improper_.insert({torsion_class, std::move(forces)});
                continue;
            }

            if (read_classes_(improper, "class", torsion_class)) {
                std::swap(torsion_class[0], torsion_class[1]);
                type_to_improper_.insert({torsion_class, std::move(forces)});
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: " << e.what() << std::endl;
        }
    }
}
