#include "./dtoa_milo.h"

#include <stdio.h>

int main() {
    char ss[32]={};
    milo::dtoa(1.2345, ss, 2);
    printf("-->%s\n", ss);s
    milo::dtoa(23.1233, ss, 21);
    printf("-->%s\n", ss);
    milo::dtoa(23.44, ss, 3);
    printf("-->%s\n", ss);
    return 0;
}