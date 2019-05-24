#ifndef INTERFACE_H
#define INTERFACE_H

#include <stddef.h>
#include "spear/spearmint_exports.hpp"

#ifdef __cplusplus
extern "C" {
#endif

SPEARMINT_EXPORT const char* spear_get_error();

/*
 * All functions return 0 upon failure and a non-zero number upon success
 */

SPEARMINT_EXPORT size_t spear_initialize_complex(const char* filename);
SPEARMINT_EXPORT size_t spear_write_complex(const char* filename);

SPEARMINT_EXPORT size_t spear_initialize_receptor(const char* filename);

SPEARMINT_EXPORT size_t spear_receptor_atom_count();
SPEARMINT_EXPORT size_t spear_receptor_atoms(float* pos);
SPEARMINT_EXPORT size_t spear_receptor_set_positions(const float* positions);
SPEARMINT_EXPORT char spear_receptor_atom_chain(size_t atom);
SPEARMINT_EXPORT size_t spear_receptor_atom_resi(size_t atom);
SPEARMINT_EXPORT char spear_receptor_atom_resn(size_t atom);
SPEARMINT_EXPORT size_t spear_receptor_atom_element(size_t atom);
SPEARMINT_EXPORT size_t spear_receptor_atom_details(
	char* cids, size_t* resi, char* resn, size_t* elements
);

SPEARMINT_EXPORT size_t spear_receptor_bond_count();
SPEARMINT_EXPORT size_t spear_receptor_bonds(size_t* bonds);
SPEARMINT_EXPORT size_t spear_receptor_bonds_in(
	size_t* atoms, size_t atoms_size, size_t* atoms_out
);

SPEARMINT_EXPORT size_t spear_initialize_ligand(const char* filename);

SPEARMINT_EXPORT size_t spear_ligand_atom_count();
SPEARMINT_EXPORT size_t spear_ligand_atoms(float* pos);
SPEARMINT_EXPORT size_t spear_ligand_set_positions(const float* positions);
SPEARMINT_EXPORT size_t spear_ligand_atom_details(size_t* elements);

SPEARMINT_EXPORT size_t spear_ligand_bond_count();
SPEARMINT_EXPORT size_t spear_ligand_bonds(size_t* bonds);
SPEARMINT_EXPORT size_t spear_ligand_bonds_in(
	size_t* atoms, size_t atoms_size, size_t* atoms_out
);
SPEARMINT_EXPORT size_t spear_ligand_neighbors(size_t atom_idx, size_t* neighbors);

SPEARMINT_EXPORT size_t spear_ligand_is_adjacent(size_t atom1, size_t atom2);
SPEARMINT_EXPORT size_t spear_ligand_add_bond(size_t atom1, size_t atom2);
SPEARMINT_EXPORT size_t spear_ligand_remove_bond(size_t atom1, size_t atom2);
SPEARMINT_EXPORT size_t spear_ligand_remove_hydrogens();
SPEARMINT_EXPORT size_t spear_ligand_add_atom(
	size_t element, float x, float y, float z
);
SPEARMINT_EXPORT size_t spear_ligand_add_atom_to(
	size_t atom, size_t element, float* x, float* y, float* z
);
SPEARMINT_EXPORT size_t spear_ligand_remove_atom(size_t atom);

SPEARMINT_EXPORT size_t spear_initialize_scoring(const char* data_dir);

SPEARMINT_EXPORT float spear_calculate_score();

#ifdef __cplusplus
}
#endif

#endif
