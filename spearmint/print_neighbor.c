#include "spearmint.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

    if (argc < 3) {
        printf("Please give 2 variables\n");
        return 1;
    }

    if (spear_initialize_ligand(argv[1]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }

    char* ptr;
    size_t idx = strtoul(argv[2], &ptr, 10);

    size_t neighbors[10];
    size_t size = spear_ligand_neighbors(idx, neighbors);

    if (size == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }

    for(size_t i = 0; i < size; ++i) {
        printf("%lu\n", neighbors[i]);
    }

    return 0;
}

