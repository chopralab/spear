#include "spearmint.h"
#include "spearmint_statics.hpp"

#include "chemfiles/Trajectory.hpp"
#include "chemfiles/Frame.hpp"

#include "spear/Molecule.hpp"
#include "spear/Constants.hpp"

size_t spear_ligand_atom_count() {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return ligand->size();
}

size_t spear_ligand_atoms(float* pos) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return get_atom_positions(*ligand, pos);
}

size_t spear_ligand_atom_details(size_t* elements) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    try {

        for (auto atom : *ligand) {
            elements[static_cast<int>(atom)] = static_cast<size_t>(atom.atomic_number());
        }

        return 1;
    } catch (std::exception& e) {
        set_error(std::string("Error creating ligand atom arrays: ") + e.what());
        return 0;
    }
}

size_t spear_ligand_bond_count() {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return ligand->topology().bonds().size();
}

size_t spear_ligand_bonds(size_t* bonds) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return get_bonds(*ligand, bonds);
}

size_t spear_ligand_neighbors(size_t atom_idx, size_t* neighbors) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    try {
        size_t current = 0;
        for (auto neighbor : (*ligand)[atom_idx].neighbors()) {
            neighbors[current++] = static_cast<size_t>(neighbor);
        }
        return current;
    } catch (std::exception& e) {
        set_error("Error creating atom arrays");
        return 0;
    }
}

size_t spear_ligand_set_positions(const float* positions) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return set_positions(*ligand, positions);
}

size_t spear_ligand_is_adjacent(size_t atom1, size_t atom2) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }
    auto& graph = ligand->graph();
    auto edge = boost::edge(atom1, atom2, graph);
    return edge.second;
}

size_t spear_ligand_add_bond(size_t atom1, size_t atom2) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }

    auto av1 = (*ligand)[atom1];
    if (av1.degree() >= av1.expected_bonds()) {
        set_error(std::string("Atom ") + std::to_string(atom1) + " is saturated.");
        return 0;
    }

    auto av2 = (*ligand)[atom2];
    if (av2.degree() >= av2.expected_bonds()) {
        set_error(std::string("Atom ") + std::to_string(atom2) + " is saturated.");
        return 0;
    }

    try {
        ligand->add_bond(atom1, atom2);
    } catch (const std::exception& e) {
        set_error(std::string("Error in adding bond ") + e.what());
        return 0;
    }
    return 1;
}

size_t spear_ligand_remove_bond(size_t atom1, size_t atom2) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }

    try {
        ligand->remove_bond(atom1, atom2);
    } catch (const std::exception& e) {
        set_error(std::string("Error in removing bond ") + e.what());
        return 0;
    }
    return 1;
}

size_t spear_ligand_remove_hydrogens() {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    try {
        ligand->remove_hydrogens();
    } catch (const std::exception& e) {
        set_error(std::string("Error in removing hydrogens ") + e.what());
        return 0;
    }
    return 1;
}

size_t spear_ligand_add_atom(size_t element, float x, float y, float z) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");

    try {
        return ligand->add_atom(static_cast<Spear::Element::Symbol>(element), {x, y, z});
    } catch (...) {
        set_error("An unknown error occured in spear_ligand_add_atom");
        return 0;
    }
}

size_t spear_ligand_add_atom_to(size_t atom, size_t element, float* x, float* y, float* z) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");

    if (x == nullptr || y == nullptr || z == nullptr) {
        set_error("X, Y, and/or Z argument(s) is invalid");
        return 0;
    }

    if (atom >= ligand->size()) {
        set_error("Atom greater than molecule size.");
        return 0;
    }

    auto av1 = (*ligand)[atom];
    if (av1.degree() >= av1.expected_bonds()) {
        set_error("Atom is saturated!");
        return 0;
    }
    auto symbol = static_cast<Spear::Element::Symbol>(element);

    try {
        auto new_idx = ligand->add_atom_to(symbol, atom);
        auto pos = ligand->positions();
        *x = pos[new_idx][0];
        *y = pos[new_idx][1];
        *z = pos[new_idx][2];
        return new_idx;
    } catch(...) {
        set_error("An unknown error occured in spear_ligand_add_atom_to");
        return 0;
    }
}
