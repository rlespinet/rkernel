#include <map>
#include <vector>
#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "mismatch.hpp"
#include "utils.hpp"

#include "debug.hpp"

template<typename hash_t>
struct kmer_data_t : kmer_count<hash_t> {
    int seq_id;
    int mismatchs;
};

template<typename hash_t>
inline bool compare(const kmer_data_t<hash_t> &s1, const kmer_data_t<hash_t> &s2) {
    return s1.seq_id < s2.seq_id;
}

template<typename hash_t>
struct depth_cache_t {
    kmer_decoding_key_t *decode_keys;
    kmer_data_t<hash_t> *kmer_data;
    int line_size;
};

template<typename seq_t, typename hash_t>
void mismatch_compute_rec(depth_cache_t<hash_t> &cache, kernel_t *kernel,
                          kmer_data_t<hash_t>* tracks, int len,
                          int k, int d, int l) {

    if (len == 0) {
        return;
    }

    if (d == k) {
        // This is a leaf !

        std::sort(tracks, tracks + len, compare<hash_t>);

        for (int i = 0; i < len - 1; i++) {
            for (int j = i + 1; j < len; j++) {

                int id1 = tracks[i].seq_id;
                int id2 = tracks[j].seq_id;

                int matches = tracks[i].count * tracks[j].count;

                kernel->data[id1 * kernel->size + id2] += matches;
                // kernel->data[id2 * kernel->size + id1] += matches;
            }
        }

        return;
    }

    kmer_data_t<hash_t> *new_tracks = cache.kmer_data + d * cache.line_size;
    kmer_decoding_key_t &key = cache.decode_keys[d];

    for (int a = 0; a < l; a++) {

        int new_len = 0;
        for (int i = 0; i < len; i++) {

            seq_t c = kmer_decode_id<seq_t, hash_t>(tracks[i], key);
            if (c == a || tracks[i].mismatchs > 0) {

                new_tracks[new_len] = tracks[i];
                new_tracks[new_len].mismatchs -= (c == a) ? 0 : 1;

                new_len++;
            }

        }

        mismatch_compute_rec<seq_t, hash_t>(cache, kernel, new_tracks, new_len, k, d + 1, l);

    }

}

template<typename seq_t, typename hash_t>
static inline bool mismatch_compute(kernel_t *kernel, kmer_data_t<hash_t>* tracks,
                                    int len, int k, int l) {


    kmer_decoding_key_t *decode_keys = new kmer_decoding_key_t[k];
    if (decode_keys == nullptr) {
        return false;
    }

    for (int i = 0; i < k; i++) {
        decode_keys[i] = kmer_decoding_key(i, k, l);
    }

    kmer_data_t<hash_t> *kmer_cache = new kmer_data_t<hash_t>[len * k];
    if (kmer_cache == nullptr) {
        delete[] decode_keys;
        return false;
    }

    depth_cache_t<hash_t> cache = {
        decode_keys,
        kmer_cache,
        len
    };

    mismatch_compute_rec<seq_t, hash_t>(cache, kernel, tracks, len, k, 0, l);

    symmetrise_lower(kernel);

    delete[] kmer_cache;
    delete[] decode_keys;

    return true;
}

kernel_t mismatch(const seq_array_t &array, int k, int m) {

    using hash_t = uint64_t;
    using seq_t = int;

    const int seq_count = array.sequences_count;
    const int seq_len = array.sequences_length;
    const int alphabet_size = array.alphabet_size;

    float *kernel_data = new float[seq_count * seq_count]();
    if (kernel_data == nullptr) {
        return invalid_kernel;
    }

    const int max_kmer_len = seq_len - k + 1;
    kmer_count<hash_t> *kmers_count = new kmer_count<hash_t>[max_kmer_len];
    if (kmers_count == nullptr) {
        delete[] kernel_data;
        return invalid_kernel;
    }

    kmer_data_t<hash_t> *kmers_data = new kmer_data_t<hash_t>[seq_count * max_kmer_len];
    if (kmers_data == nullptr) {
        delete[] kmers_count;
        delete[] kernel_data;
        return invalid_kernel;
    }

    // Filling in the structure
    int kmers_len = 0;
    for (int i = 0; i < seq_count; i++) {
        int len = kmer_encoding_count(kmers_count, array.data + i * seq_len,
                                      seq_len, k, alphabet_size);

        for (int j = 0; j < len; j++) {
            (kmer_count<hash_t>&) kmers_data[kmers_len] = kmers_count[j];

            kmers_data[kmers_len].seq_id = i;
            kmers_data[kmers_len].mismatchs = m;

            kmers_len++;
        }

    }

    kernel_t kernel = {
        kernel_data,
        seq_count
    };

    bool ok = mismatch_compute<seq_t, hash_t>(&kernel, kmers_data, kmers_len, k, alphabet_size);
    if (!ok) {
        delete[] kmers_data;
        delete[] kmers_count;
        delete[] kernel_data;
        return invalid_kernel;
    }

    delete[] kmers_data;
    delete[] kmers_count;

    return kernel;
}
