#pragma once

#include "common.hpp"

template<typename letter, typename dtype>
sq_matrix<dtype> spectrum(const vector2D<letter> &sequences,
                          int sequences_len, int alphabet_size, int k);

template<typename letter, typename dtype>
matrix<dtype> spectrum(const vector2D<letter> &sequences1,
                       const vector2D<letter> &sequences2,
                       int sequences_len, int alphabet_size, int k);
