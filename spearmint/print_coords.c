#include "spearmint.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    if (argc < 2) {
        printf("Please give 1 variable\n");
        return 1;
    }

    if (spear_initialize_ligand(argv[1]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }

    size_t ligand_size = spear_ligand_atom_count();

    float* pos = (float*)malloc(sizeof(float)*ligand_size*3);
    spear_ligand_atoms(pos);

    size_t* atn = (size_t*)malloc(sizeof(size_t)*ligand_size);
    spear_ligand_atom_details(atn);

    printf("%s\n", argv[1]);
    printf("%lu\n", ligand_size);
    for(size_t i = 0; i < ligand_size; ++i) {
        printf("%lu\t%f\t%f\t%f\n", atn[i], pos[i*3+0], pos[i*3+1], pos[i*3+2]);
    }

    free(pos);
    free(atn);

    return 0;
}

