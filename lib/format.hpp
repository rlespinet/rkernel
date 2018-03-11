#pragma once

#include "common.hpp"

template<typename letter>
struct sequence_array {
    int alphabet_size;
    int sequences_len;
    vector2D<letter> sequences;
};

template<typename letter, typename dtype>
sequence_array<letter> to_sequence_array(dtype *data, int rows, int cols);
