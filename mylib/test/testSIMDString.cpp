#include "../mystl/string.h"
#include <stdio.h>

int main(){
    SIMDString<> ss="@@@!!!", s2="154679",s3="6.94648";
    
    printf("%d,%lf,%s,%s\n",stoi(s2),stof(s3),
    to_string(9876543).c_str(),to_string(9998.98765).c_str());
}