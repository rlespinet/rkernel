#include <map>
#include <vector>
#include <cstring>

#include "spectrum.hpp"
#include "utils.hpp"

#include "debug.hpp"

template<class K, class V>
using map = std::map<K, V>;

struct mismatch_struct_t {
    const int *ptr;
    int seq_id;
};

static void compute_kernel_rec(kernel_t &kernel, mismatch_struct_t* tracks,
                               int *mismatchs, int len, int k, int d, int l) {

    // static std::vector<int> debug;

    if (len == 0) {
        return;
    }

    if (d == k) {
        // This is a leaf !
        for (int i = 0; i < len - 1; i++) {
            for (int j = i + 1; j < len; j++) {

                int id1 = tracks[i].seq_id;
                int id2 = tracks[j].seq_id;

                kernel.data[id1 * kernel.size + id2] += 1;
                kernel.data[id2 * kernel.size + id1] += 1;
            }
        }

        return;
    }


    int *old_mismatchs = new int[len];
    mismatch_struct_t *old_tracks = new mismatch_struct_t[len];

    for (int a = 0; a < l; a++) {

        // debug.push_back(a);
        // for (int i = 0; i < debug.size(); i++) {
        //     std::cout << debug[i] << " ";
        // }

        std::memcpy(old_mismatchs, mismatchs, sizeof(int) * len);
        std::memcpy(old_tracks, tracks, sizeof(mismatch_struct_t) * len);

        int i = 0, j = len;
        while (i < j) {

            if (tracks[i].ptr[d] == a || mismatchs[i] > 0) {

                if (tracks[i].ptr[d] != a) {
                    mismatchs[i]--;
                }

                i++;
            } else {
                j--;
                std::swap(tracks[i], tracks[j]);
                std::swap(mismatchs[i], mismatchs[j]);
            }

        }

        // std::cout << ": " << j << std::endl;

        compute_kernel_rec(kernel, tracks, mismatchs, j, k, d + 1, l);

        std::memcpy(mismatchs, old_mismatchs, sizeof(int) * len);
        std::memcpy(tracks, old_tracks, sizeof(mismatch_struct_t) * len);

        // debug.pop_back();

    }

}


kernel_t mismatch(const seq_array_t &array, int k, int m) {

    const int *data = array.data;
    const int rows = array.sequences_count;
    const int cols = array.sequences_length;
    const int l = array.alphabet_size;

    float *kernel_data = new float[rows * rows]();
    if (kernel_data == nullptr) {
        return invalid_kernel;
    }

    int tracks_len = rows * (cols - k + 1);
    mismatch_struct_t *tracks = new mismatch_struct_t[tracks_len];
    if (tracks == nullptr) {
        return invalid_kernel;
    }

    for (int i = 0; i < rows; i++) {

        const int* data_i = data + cols * i;

        for (int j = 0; j < cols - k + 1; j++) {

            tracks[i * (cols - k + 1) + j].ptr = data_i + j;
            tracks[i * (cols - k + 1) + j].seq_id = i;

        }

    }

    int *mismatchs = new int[tracks_len];
    if (mismatchs == nullptr) {
        return invalid_kernel;
    }


    for (int i = 0; i < tracks_len; i++) {
        mismatchs[i] = m;
    }

    kernel_t kernel = {
        kernel_data,
        rows
    };

    compute_kernel_rec(kernel, tracks, mismatchs, tracks_len, k, 0, l);

    delete[] tracks;
    delete[] mismatchs;

    return kernel;
}
