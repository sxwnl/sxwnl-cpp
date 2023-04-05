#ifndef MYTINYSTL_ASTRING_H_
#define MYTINYSTL_ASTRING_H_

#include "string/SIMDString.h"

namespace mystl {
    // template<size_t INTERNAL_SIZE = 64, class Allocator = ::std::allocator<char>>
    // using string = SIMDString<INTERNAL_SIZE, Allocator>;
    using string = SIMDString<>;
} // namespace mystl
#endif // !MYTINYSTL_ASTRING_H_