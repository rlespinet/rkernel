#pragma once

#include "common.hpp"

template<typename dtype>
sq_matrix<dtype> local_alignment(const vector2D<ltype> &sequences,
                                 int sequences_len, int alphabet_size,
                                 dtype substitution, dtype gap_open,
                                 dtype gap_extension);
