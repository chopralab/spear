#include "spearmint.h"
#include <stdio.h>

int main(int argc, char** argv) {

    if (argc < 3) {
        printf("Please give 2 variables\n");
        return 1;
    }

    if (spear_initialize_complex(argv[1]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }
    if (spear_write_complex(argv[2]) == 0) {
        printf("%s\n", spear_get_error());
        return 1;
    }

    return 0;
}
