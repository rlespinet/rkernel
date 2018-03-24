#pragma once

#include "common.hpp"

template<typename dtype>
sq_matrix<dtype> needleman_wunsch_affine_gap(const vector2D<ltype> &sequences,
                                             int sequences_len, int alphabet_size,
                                             dtype similarity, dtype substitution,
                                             dtype gap_open, dtype gap_extension);

template<typename dtype>
sq_matrix<dtype> smith_waterman_affine_gap(const vector2D<ltype> &sequences,
                                           int sequences_len, int alphabet_size,
                                           dtype similarity, dtype substitution,
                                           dtype gap_open, dtype gap_extension);
