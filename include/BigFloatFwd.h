#pragma once

// This file is for if only a forward declaration of BigFloat is needed,
// e.g. if a file only has a pointer or reference to the type and it is
// not actually used in the file.

#include <stddef.h>

namespace big_float {

template<size_t THE_N>
struct BigFloat;

}
