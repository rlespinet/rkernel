#pragma once

#include "common.hpp"

struct sequence_array {
    int alphabet_size;
    int sequences_len;
    vector2D<ltype> sequences;
};

template<typename dtype>
sequence_array to_sequence_array(dtype *data, int rows, int cols);
