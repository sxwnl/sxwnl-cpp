#include "./dtoa_milo.h"


int main() {
    char ss[128]={};
    milo::dtoa(65448566564654, ss);
    printf("-->%s\n", ss);
    return 0;
}