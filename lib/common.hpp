#pragma once

#include <cstdint>

struct seq_array_t {
    int *data;
    int sequences_count;
    int sequences_length;
    int alphabet_size;
};

const seq_array_t invalid_seq_array = {nullptr, 0, 0, 0};

inline bool is_valid(const seq_array_t &seq_array) {
    return seq_array.data != nullptr;
}

inline bool is_invalid(const seq_array_t &seq_array) {
    return !is_valid(seq_array);
}

struct kernel_t {
    float *data;
    int size;
};

inline bool is_valid(const kernel_t &kernel) {
    return kernel.data != nullptr;
}

inline bool is_invalid(const kernel_t &kernel) {
    return !is_valid(kernel);
}

const kernel_t invalid_kernel = {nullptr, 0};
