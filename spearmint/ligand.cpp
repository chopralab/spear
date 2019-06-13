#include "spearmint.h"
#include "spearmint_statics.hpp"

#include "chemfiles/Trajectory.hpp"
#include "chemfiles/Frame.hpp"

#include "spear/Molecule.hpp"
#include "spear/Constants.hpp"

uint64_t spear_ligand_atom_count() {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return ligand->size();
}

uint64_t spear_ligand_atoms(float* pos) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return get_atom_positions(*ligand, pos);
}

uint64_t spear_ligand_atom_details(uint64_t* elements) {
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

uint64_t spear_ligand_bond_count() {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return ligand->topology().bonds().size();
}

uint64_t spear_ligand_bonds(uint64_t* bonds) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return get_bonds(*ligand, bonds);
}

uint64_t spear_ligand_bonds_in(uint64_t* atoms, uint64_t atoms_size, uint64_t* atoms_out) {
	CHECK_MOLECULE(ligand, "You must run initialize_receptor first");
	std::set<size_t> atoms_to_check(atoms, atoms + atoms_size);

	auto res = ligand->get_bonds_in(atoms_to_check);

	for (size_t i = 0; i < res.size(); ++i) {
		atoms_out[3 * i + 0] = res[i].source();
		atoms_out[3 * i + 1] = res[i].target();
		atoms_out[3 * i + 2] = res[i].order();
	}

	return res.size();
}

uint64_t spear_ligand_neighbors(uint64_t atom_idx, uint64_t* neighbors) {
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

uint64_t spear_ligand_set_positions(const float* positions) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    return set_positions(*ligand, positions);
}

uint64_t spear_ligand_is_adjacent(uint64_t atom1, uint64_t atom2) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }
    auto& graph = ligand->graph();
    auto edge = boost::edge(atom1, atom2, graph);
    return edge.second;
}

uint64_t spear_ligand_add_bond(uint64_t atom1, uint64_t atom2) {
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
    } catch (...) {
        set_error("Unknown error in spear_ligand_add_bond");
        return 0;
    }
    return 1;
}

uint64_t spear_ligand_remove_bond(uint64_t atom1, uint64_t atom2) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    if (atom1 > ligand->size() || atom2 > ligand->size()) {
        return 0;
    }

    try {
        ligand->remove_bond(atom1, atom2);
    } catch (const std::exception& e) {
        set_error(std::string("Error in removing bond ") + e.what());
        return 0;
    } catch (...) {
        set_error("Unknown spear_ligand_remove_bond");
        return 0;
    }
    return 1;
}

uint64_t spear_ligand_remove_hydrogens() {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");
    try {
        ligand->remove_hydrogens();
    } catch (const std::exception& e) {
        set_error(std::string("Error in removing hydrogens ") + e.what());
        return 0;
    } catch (...) {
        set_error("Unknown error in spear_ligand_remove_hydrogens");
        return 0;
    }
    return 1;
}

uint64_t spear_ligand_add_atom(uint64_t element, float x, float y, float z) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");

    try {
        return ligand->add_atom(static_cast<Spear::Element::Symbol>(element), {x, y, z});
    } catch (...) {
        set_error("An unknown error occured in spear_ligand_add_atom");
        return 0;
    }
}

uint64_t spear_ligand_add_atom_to(uint64_t atom, uint64_t element, float* x, float* y, float* z) {
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

uint64_t spear_ligand_swap_atoms(uint64_t idx1, uint64_t idx2) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");

    try {
        ligand->swap_atoms(idx1, idx2);
        return 1;
    } catch(const std::exception& error) {
        set_error(error.what());
    } catch(...) {
        set_error("An unknown error occured in spear_ligand_swap_atoms");
    }

    return 0;
}

uint64_t spear_ligand_remove_atom(uint64_t atom) {
    CHECK_MOLECULE(ligand, "You must run initialize_ligand first");

    if (atom >= ligand->size()) {
        set_error("Atom not in range for removal.");
        return 0;
    }

    try {
        ligand->remove_atom(atom);
        return 1;
    } catch(...) {
        set_error("Unknown error in atom removal.");
        return 0;
    }
}