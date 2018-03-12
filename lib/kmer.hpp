#pragma once

#include "common.hpp"
#include "utils.hpp"

#include <map>
#include <vector>

template<class K, class V>
using tree_map = std::map<K, V>;

// Type in which we encode kmers
using e_type = uint64_t;

struct kmer {
    e_type encoding;
    int count;
};

inline int compare(const kmer &s1, const kmer &s2) {
    return (int) (s1.encoding - s2.encoding);
}

template<typename letter>
vector1D<kmer> kmer_encode(const view1D<letter> &sequence, int k, int sequence_len, int alphabet_size) {

    const e_type shift_mod = ipow(alphabet_size, k - 1);

    e_type encoding = 0;
    for (int j = 0; j < k - 1; j++) {
        encoding = encoding * alphabet_size + static_cast<e_type>(sequence[j]);
    }

    tree_map<e_type, int> map;

    for (int j = 0; j < sequence_len - k + 1; j++) {
        encoding = encoding * alphabet_size + static_cast<e_type>(sequence[j + k - 1]);

        map[encoding]++;

        encoding = encoding % shift_mod;
    }

    int count = 0;
    vector1D<kmer> kmers(sequence_len - k + 1);
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        kmers[count++] = {it->first, it->second};
    }
    kmers.resize(count);

    return kmers;
}

template<typename letter>
vector2D<kmer> kmer_encode_all(const vector2D<letter> &sequences, int k,
                               int sequences_len, int alphabet_size) {

    vector2D<kmer> all_kmers(sequences.size(), sequences_len);

    for (int i = 0; i < sequences.size(); i++) {
        all_kmers[i] = kmer_encode(sequences[i], k, sequences_len, alphabet_size);
    }

    return all_kmers;
}
