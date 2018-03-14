#pragma once

#include "common.hpp"

template<typename dtype>
sq_matrix<dtype> spectrum(const vector2D<ltype> &sequences,
                          int sequences_len, int alphabet_size, int k);

template<typename dtype>
matrix<dtype> spectrum(const vector2D<ltype> &sequences1,
                       const vector2D<ltype> &sequences2,
                       int sequences_len, int alphabet_size, int k);
