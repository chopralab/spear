#include "spearmint.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    if (argc < 4) {
        printf("Please give 3 variables\n");
        return 1;
    }

    if (spear_receptor_atom_count() == 0) {
        printf("(Ignore me: %s)\n", spear_get_error());
    } else {
        printf("Something really bad happened\n");
    }

    if (spear_ligand_atom_count() == 0) {
        printf("(Ignore me: %s)\n", spear_get_error());
    } else {
        printf("Something really bad happened\n");
    }

    if (spear_initialize_receptor(argv[1]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }
    if (spear_initialize_ligand(argv[2]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }
    if (spear_initialize_scoring(argv[3]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }

    float score = spear_calculate_score();

    printf("score is %f\n", score);

    size_t lig_atom_count = spear_ligand_atom_count();

    printf("There's %lu ligand atoms\n", lig_atom_count);
    printf("There's %lu receptor atoms\n", spear_receptor_atom_count());

    float x, y, z;

    spear_ligand_add_atom_to(12, 6, &x, &y, &z);
    score = spear_calculate_score();
    printf("After methyl addition at (%f,%f,%f), score is %f\n", score, x, y, z);

    float* pos = (float*)malloc(sizeof(float)*spear_ligand_atom_count()*3);
    spear_ligand_atoms(pos);

    for (size_t i = 0; i < spear_ligand_atom_count(); ++i) {
        // Change the X position of atom 'i' by 5 angstroms
        pos[ i * 3 + 0 ] = pos[ i * 3 + 0 ] + 5.0f;
                        
        // Change the Y position of atom 'i' by 4 angstroms
        pos[ i * 3 + 1 ] = pos[ i * 3 + 1 ] + 4.0f;
                        
        // Change the Z position of atom 'i' by 3 angstroms
        pos[ i * 3 + 2 ] = pos[ i * 3 + 2 ] + 3.0f;
    }

    if (spear_ligand_set_positions(pos) == 0) {
        printf("%s\n", spear_get_error());
    }

    score = spear_calculate_score();
    // Nothing else needs to be done to recalculate the score.
    printf("After moving, score is %f\n", score);

    free(pos);
    spear_write_complex("out.pdb");

    return 0;
}
