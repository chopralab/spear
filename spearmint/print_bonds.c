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

    size_t bond_size = spear_ligand_bond_count();

    size_t* pos = (size_t*)malloc(sizeof(size_t)*bond_size*3);
    spear_ligand_bonds(pos);

    for(size_t i = 0; i < bond_size; ++i) {
        printf("%lu\t%lu\t%lu\n", pos[i*3+0], pos[i*3+1], pos[i*3+2]);
    }

    free(pos);

    return 0;
}

