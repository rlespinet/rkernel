#pragma once

#include "common.hpp"
#include "utils.hpp"

#include <map>

template<class K, class V>
using tree_map = std::map<K, V>;

template<typename hash_t>
struct kmer_count {
    hash_t hash;
    int count;
};

template<typename seq_t, typename hash_t>
inline int kmer_encoding_count(kmer_count<hash_t>* kmers, const seq_t *sequence, int seq_len, int k, int l) {

    const hash_t shift_mod = ipow(l, k - 1);

    hash_t hash = 0;
    for (int j = 0; j < k - 1; j++) {
        hash = hash * l + static_cast<hash_t>(sequence[j]);
    }

    tree_map<hash_t, int> map;

    for (int j = 0; j < seq_len - k + 1; j++) {
        hash = hash * l + static_cast<hash_t>(sequence[j + k - 1]);

        map[hash]++;

        hash = hash % shift_mod;
    }

    int size = 0;
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        kmers[size].hash = it->first;
        kmers[size].count = it->second;
        size++;
    }

    return size;

}

struct kmer_decoding_key_t {
    int div;
    int l;
};

inline kmer_decoding_key_t kmer_decoding_key(int i, int k, int l) {
    kmer_decoding_key_t key = {
        ipow(l, k - 1 - i), l
    };
    return key;
}

template<typename seq_t, typename hash_t>
inline seq_t kmer_decode_id(const kmer_count<hash_t> &kmers, const kmer_decoding_key_t &key) {
    return (kmers.hash / key.div) % key.l;
}
