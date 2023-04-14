//#include "../mystl/my_string.h"

#include <cstdio>
#include <string>
#include <limits>
#include <cstdlib>
int main() {
    char ss[32] = {};
    printf("%s\n", ss);
    printf("%d\n", std::numeric_limits<int>::digits10 + 3);
    return 0;
}