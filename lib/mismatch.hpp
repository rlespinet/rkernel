#pragma once

#include "common.hpp"

template<typename letter, typename dtype>
sq_matrix<dtype> mismatch(const vector2D<letter> &sequences,
                          int sequences_len, int alphabet_size,
                          int k, int m);
