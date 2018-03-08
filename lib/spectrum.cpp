#include <limits>
#include <map>

#include "spectrum.hpp"
#include "utils.hpp"

#include "debug.hpp"

template<class K, class V>
using map = std::map<K, V>;

#if 1

kernel_t spectrum(const seq_array_t &array, int k) {

    const int *data = array.data;
    const int rows = array.sequences_count;
    const int cols = array.sequences_length;
    const int l = array.alphabet_size;

    if (k > ilog(std::numeric_limits<int>::max(), l)) {
        // TODO(RL) with the method used, We are in trouble in this
        // case because we cannot encode sequences of size k with
        // alphabet_size l in one int... We should warn the user about
        // what happened though...
        return invalid_kernel;
    }


    float *kernel_data = new float[rows * rows];
    if (kernel_data == nullptr) {
        return invalid_kernel;
    }

    map<int, int> *count_kmers = new map<int, int>[rows];

    int shift_mod = ipow(l, k - 1);

    for (int i = 0; i < rows; i++) {

        const int* data_i = data + cols * i;

        map<int, int> &map_i = count_kmers[i];

        int id = 0;
        for (int j = 0; j < k - 1; j++) {
            id = id * l + data_i[j];
        }

        for (int j = k - 1; j < cols; j++) {
            id = id * l + data_i[j];

            map_i[id]++;

            id = id % shift_mod;

        }
    }

    int max_kmers = (cols - k + 1);

    std::pair<int, int> *seq_kmers = new std::pair<int, int>[rows * max_kmers]();
    if (seq_kmers == nullptr) {
        return invalid_kernel;
    }

    for (int i = 0; i < rows; i++) {

        int j = 0;
        for (auto it = count_kmers[i].cbegin();
             it != count_kmers[i].cend(); ++it) {
            seq_kmers[i * max_kmers + j++] = *it;
        }

    }

    delete[] count_kmers;

    // TODO(RL) Do something smarter than omp parallel for...
// #pragma omp parallel for
    for (int i = 0; i < rows; i++) {

        const std::pair<int, int> *kmers_i = seq_kmers + i * max_kmers;

        for (int j = i; j < rows; j++) {

            const std::pair<int, int> *kmers_j = seq_kmers + j * max_kmers;

            int matches = 0;
            for (int ki = 0, kj = 0; ki < max_kmers && kj < max_kmers;) {
                if (kmers_i[ki].first < kmers_j[kj].first) {
                    ki++;
                } else if (kmers_i[ki].first > kmers_j[kj].first) {
                    kj++;
                } else {
                    matches += kmers_i[ki].second * kmers_j[kj].second;
                    ki++;
                    kj++;
                }
            }

            kernel_data[i * rows + j] = matches;
            kernel_data[j * rows + i] = matches;

        }

    }

    delete[] seq_kmers;

    kernel_t kernel = {
        kernel_data,
        rows
    };

    return kernel;
}

#else

kernel_t spectrum(const seq_array_t &array, int k) {

    const int *data = array.data;
    const int rows = array.sequences_count;
    const int cols = array.sequences_length;
    const int l = array.alphabet_size;

    float *kernel_data = new float[rows * rows];
    if (kernel_data == nullptr) {
        return invalid_kernel;
    }

    map<int, int> *count_kmers = new map<int, int>[rows];

    int shift_mod = ipow(l, k - 1);

    for (int i = 0; i < rows; i++) {

        const int* data_i = data + cols * i;

        map<int, int> &map_i = count_kmers[i];

        int id = 0;
        for (int j = 0; j < k - 1; j++) {
            id = id * l + data_i[j];
        }

        for (int j = k - 1; j < cols; j++) {
            id = id * l + data_i[j];

            map_i[id]++;

            id = id % shift_mod;

        }
    }


    for (int i = 0; i < rows; i++) {

        const map<int, int> &map_i = count_kmers[i];

        for (int j = i; j < rows; j++) {

            int matches = 0;

            const map<int, int> &map_j = count_kmers[j];

            for (auto it1 = map_i.cbegin(); it1 != map_i.cend(); ++it1) {

                auto it2 = map_j.find(it1->first);
                if (it2 != map_j.end()) {
                    matches += it1->second * it2->second;
                }

            }

            kernel_data[i * rows + j] = matches;
            kernel_data[j * rows + i] = matches;

        }

    }

    delete[] count_kmers;

    kernel_t kernel = {
        kernel_data,
        rows
    };

    return kernel;



}

#endif
