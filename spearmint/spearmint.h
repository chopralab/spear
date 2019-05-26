#ifndef INTERFACE_H
#define INTERFACE_H

#include <stddef.h>
#include <stdint.h>
#include "spear/spearmint_exports.hpp"

#ifdef __cplusplus
extern "C" {
#endif

SPEARMINT_EXPORT const char* spear_get_error();

/*
 * All functions return 0 upon failure and a non-zero number upon success
 */

SPEARMINT_EXPORT uint64_t spear_initialize_complex(const char* filename);
SPEARMINT_EXPORT uint64_t spear_write_complex(const char* filename);

SPEARMINT_EXPORT uint64_t spear_initialize_receptor(const char* filename);

SPEARMINT_EXPORT uint64_t spear_receptor_atom_count();
SPEARMINT_EXPORT uint64_t spear_receptor_atoms(float* pos);
SPEARMINT_EXPORT uint64_t spear_receptor_set_positions(const float* positions);
SPEARMINT_EXPORT uint8_t spear_receptor_atom_chain(uint64_t atom);
SPEARMINT_EXPORT uint64_t spear_receptor_atom_resi(uint64_t atom);
SPEARMINT_EXPORT uint8_t spear_receptor_atom_resn(uint64_t atom);
SPEARMINT_EXPORT uint64_t spear_receptor_atom_element(uint64_t atom);
SPEARMINT_EXPORT uint64_t spear_receptor_atom_details(
	uint8_t* cids, uint64_t* resi, uint8_t* resn, uint64_t* elements
);

SPEARMINT_EXPORT uint64_t spear_receptor_bond_count();
SPEARMINT_EXPORT uint64_t spear_receptor_bonds(uint64_t* bonds);
SPEARMINT_EXPORT uint64_t spear_receptor_bonds_in(
	uint64_t* atoms, uint64_t atoms_size, uint64_t* atoms_out
);

SPEARMINT_EXPORT uint64_t spear_initialize_ligand(const char* filename);

SPEARMINT_EXPORT uint64_t spear_ligand_atom_count();
SPEARMINT_EXPORT uint64_t spear_ligand_atoms(float* pos);
SPEARMINT_EXPORT uint64_t spear_ligand_set_positions(const float* positions);
SPEARMINT_EXPORT uint64_t spear_ligand_atom_details(uint64_t* elements);

SPEARMINT_EXPORT uint64_t spear_ligand_bond_count();
SPEARMINT_EXPORT uint64_t spear_ligand_bonds(uint64_t* bonds);
SPEARMINT_EXPORT uint64_t spear_ligand_bonds_in(
	uint64_t* atoms, uint64_t atoms_size, uint64_t* atoms_out
);
SPEARMINT_EXPORT uint64_t spear_ligand_neighbors(uint64_t atom_idx, uint64_t* neighbors);

SPEARMINT_EXPORT uint64_t spear_ligand_is_adjacent(uint64_t atom1, uint64_t atom2);
SPEARMINT_EXPORT uint64_t spear_ligand_add_bond(uint64_t atom1, uint64_t atom2);
SPEARMINT_EXPORT uint64_t spear_ligand_remove_bond(uint64_t atom1, uint64_t atom2);
SPEARMINT_EXPORT uint64_t spear_ligand_remove_hydrogens();
SPEARMINT_EXPORT uint64_t spear_ligand_add_atom(
	uint64_t element, float x, float y, float z
);
SPEARMINT_EXPORT uint64_t spear_ligand_add_atom_to(
	uint64_t atom, uint64_t element, float* x, float* y, float* z
);
SPEARMINT_EXPORT uint64_t spear_ligand_remove_atom(uint64_t atom);

SPEARMINT_EXPORT uint64_t spear_initialize_scoring(const char* data_dir);

SPEARMINT_EXPORT float spear_calculate_score();

#ifdef __cplusplus
}
#endif

#endif
