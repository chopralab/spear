#include "spearmint.h"
#include <stdio.h>

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

    return 0;
}
