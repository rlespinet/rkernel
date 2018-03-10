#include <limits>
#include <map>

#include "spectrum.hpp"
#include "utils.hpp"
#include "kmer.hpp"

#include "debug.hpp"

template<class K, class V>
using map = std::map<K, V>;

kernel_t spectrum(const seq_array_t &array, int k) {

    const int *data = array.data;
    const int rows = array.sequences_count;
    const int cols = array.sequences_length;
    const int l = array.alphabet_size;

    float *kernel_data = new float[rows * rows];
    if (kernel_data == nullptr) {
        return invalid_kernel;
    }

    int *kmers_len = new int[rows];
    if (kmers_len == nullptr) {
        return invalid_kernel;
    }

    const int max_kmer_len = cols - k + 1;
    kmer_count<uint64_t> *kmers = new kmer_count<uint64_t>[rows * max_kmer_len];
    if (kmers == nullptr) {
        return invalid_kernel;
    }

    for (int i = 0; i < rows; i++) {
        kmers_len[i] = kmer_encoding_count(kmers + i * max_kmer_len, &data[i * cols], cols, k, l);
    }

    // TODO(RL) Do something smarter than omp parallel for...
#pragma omp parallel for
    for (int i = 0; i < rows; i++) {

        const auto *kmers_i = kmers + i * max_kmer_len;
        const int kmers_len_i = kmers_len[i];

        for (int j = i; j < rows; j++) {

            const auto *kmers_j = kmers + j * max_kmer_len;
            const int kmers_len_j = kmers_len[j];

            int matches = 0;
            for (int ki = 0, kj = 0; ki < kmers_len_i && kj < kmers_len_j;) {
                if (kmers_i[ki].hash < kmers_j[kj].hash) {
                    ki++;
                } else if (kmers_i[ki].hash > kmers_j[kj].hash) {
                    kj++;
                } else {
                    matches += kmers_i[ki].count * kmers_j[kj].count;
                    ki++;
                    kj++;
                }
            }

            kernel_data[i * rows + j] = matches;
            kernel_data[j * rows + i] = matches;

        }

    }

    delete[] kmers;
    delete[] kmers_len;

    kernel_t kernel = {
        kernel_data,
        rows
    };

    return kernel;
}
