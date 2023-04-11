#include "../mystl/my_string/dtoa.h"

#include <stdio.h>

int main() {
    char ss[32] = {};
    dtoa_milo2(987.6666666, ss, 16, 0);
    printf("%s\n", ss);
    return 0;
}