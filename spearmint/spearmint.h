#ifndef INTERFACE_H
#define INTERFACE_H

#include <stddef.h>
#include "spear/exports.hpp"

#ifdef __cplusplus
extern "C" {
#endif

SPEAR_EXPORT const char* spear_get_error();

/*
 * All functions return 0 upon failure and a non-zero number upon success
 */

SPEAR_EXPORT size_t spear_initialize_complex(const char* filename);
SPEAR_EXPORT size_t spear_write_complex(const char* filename);

SPEAR_EXPORT size_t spear_initialize_receptor(const char* filename);

SPEAR_EXPORT size_t spear_receptor_atom_count();
SPEAR_EXPORT size_t spear_receptor_atoms(float* pos);
SPEAR_EXPORT size_t spear_receptor_set_positions(const float* positions);
SPEAR_EXPORT size_t spear_receptor_atom_details(char* cids, size_t* resi,
                                                char* resn, size_t* elements);

SPEAR_EXPORT size_t spear_receptor_bond_count();
SPEAR_EXPORT size_t spear_receptor_bonds(size_t* bonds);

SPEAR_EXPORT size_t spear_initialize_ligand(const char* filename);

SPEAR_EXPORT size_t spear_ligand_atom_count();
SPEAR_EXPORT size_t spear_ligand_atoms(float* pos);
SPEAR_EXPORT size_t spear_ligand_set_positions(const float* positions);
SPEAR_EXPORT size_t spear_ligand_atom_details(size_t* elements);

SPEAR_EXPORT size_t spear_ligand_bond_count();
SPEAR_EXPORT size_t spear_ligand_bonds(size_t* bonds);
SPEAR_EXPORT size_t spear_ligand_neighbors(size_t atom_idx, size_t* neighbors);

SPEAR_EXPORT size_t spear_ligand_is_adjacent(size_t atom1, size_t atom2);
SPEAR_EXPORT size_t spear_ligand_add_bond(size_t atom1, size_t atom2);
SPEAR_EXPORT size_t spear_ligand_remove_bond(size_t atom1, size_t atom2);
SPEAR_EXPORT size_t spear_ligand_remove_hydrogens();
SPEAR_EXPORT size_t spear_ligand_add_atom(size_t element, float x, float y, float z);
SPEAR_EXPORT size_t spear_ligand_add_atom_to(size_t atom, size_t element, float* x, float* y, float* z);

SPEAR_EXPORT size_t spear_initialize_scoring(const char* data_dir);

SPEAR_EXPORT float spear_calculate_score();

#ifdef __cplusplus
}
#endif

#endif
