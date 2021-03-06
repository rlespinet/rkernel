#pragma once

#include "common.hpp"

template<typename dtype>
sq_matrix<dtype> mismatch(const vector2D<ltype> &sequences,
                          int sequences_len, int alphabet_size,
                          int k, int m);
