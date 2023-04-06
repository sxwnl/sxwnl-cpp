#include "./dtoa_milo.h"

#include <stdio.h>

int main() {
    char ss[32]={};
    
    dtoa_milo2(0.333333333333, ss, 5, 1);
    printf("-->%s\n", ss);
    dtoa_milo(0.123123, ss);
    printf("-->%s\n", ss);
    double v = 123.123456789;
    int len, K;
    DtoaMilo::Grisu2(v, ss, &len, &K);
    DtoaMilo::Prettify(ss, len, K, 5, 0);
    printf("-->%s\n", ss);

    return 0;
}