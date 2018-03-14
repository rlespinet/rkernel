#pragma once

#include "common.hpp"
#include "utils.hpp"

#include <unordered_map>
#include <map>

template<class K, class V>
using tree_map = std::map<K, V>;

template<class K, class V>
using hash_map = std::unordered_map<K, V>;

struct kmer {
    const ltype *data;
    int k;

    kmer(const ltype *seq = nullptr, int size = 0)
        : data(seq)
        , k(size) {
    }

    kmer(const kmer& oth)
        : data(oth.data)
        , k(oth.k) {
    }

    kmer(kmer&& oth) {
        swap(*this, oth);
    }

    ~kmer() {
    }

    kmer& operator=(kmer oth) {
        swap(*this, oth);
        return *this;
    }

    void swap(kmer& s1, kmer& s2) {
        using std::swap;

        swap(s1.data, s2.data);
        swap(s1.k, s2.k);
    }

    const ltype& operator[](int i) const {
        return data[i];
    }

    int size() const {
        return k;
    }

};

inline int compare(const kmer &s1, const kmer &s2) {
    for (int i = 0; ; i++) {
        if (i >= std::min(s1.size(), s2.size())) {
            return s1.size() - s2.size();
        }
        if (s1[i] != s2[i]) {
            return s1[i] - s2[i];
        }
    }
}

inline bool operator<(const kmer &s1, const kmer &s2) {
    return compare(s1, s2) < 0;
}

inline bool operator>(const kmer &s1, const kmer &s2) {
    return compare(s1, s2) > 0;
}

inline bool operator==(const kmer &s1, const kmer &s2) {
    if (s1.size() != s2.size()) {
        return false;
    }
    for (int i = 0; i < s1.size(); i++) {
        if (s1[i] != s2[i]) {
            return false;
        }
    }
    return true;
}

inline bool operator>=(const kmer &s1, const kmer &s2) {
    return !(s1 < s2);
}

inline bool operator<=(const kmer &s1, const kmer &s2) {
    return !(s1 > s2);
}

inline bool operator!=(const kmer &s1, const kmer &s2) {
    return !(s1 == s2);
}

namespace std {

    template <>
    struct hash< kmer >
    {
        std::size_t operator()(const kmer& s) const {
            using std::hash;
            using std::size_t;

            size_t h = 0;
            for (int i = 0; i < s.size(); i++) {
                h = h * 709 + s[i];
            }
            return h;
        }
    };

}

struct kmer_count {
    kmer data;
    int count;

    kmer_count(const kmer &s, int c)
        : data(s)
        , count(c) {
    }

    kmer_count(const ltype *seq = nullptr, int size = 0, int c = 0)
        : data(seq, size)
        , count(c) {
    }

    ~kmer_count() {
    }

};

inline vector1D< kmer_count > count_kmers(const vector1D<ltype> &seq, int k) {

    tree_map<kmer, int> map;

    for (int i = 0; i < seq.size() - k + 1; i++) {
        kmer s(seq.data() + i, k);
        map[s]++;
    }

    int id = 0;
    vector1D< kmer_count > result(map.size());
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        result[id++] = kmer_count(it->first, it->second);
    }

    return result;
}

inline vector2D< kmer_count > count_all_kmers(const vector2D<ltype> &seqs, int k) {

    vector2D< kmer_count > result(seqs.size(), 0);

    for (int i = 0; i < seqs.size(); i++) {
        result[i] = count_kmers(seqs[i], k);
    }

    return result;
}
