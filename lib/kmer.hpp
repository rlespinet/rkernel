#pragma once

#include "common.hpp"
#include "utils.hpp"

#include <unordered_map>
#include <map>

template<class K, class V>
using tree_map = std::map<K, V>;

template<class K, class V>
using hash_map = std::unordered_map<K, V>;

template<typename letter>
struct kmer {
    const letter *data;
    int k;

    kmer(const letter *seq = nullptr, int size = 0)
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

    const letter& operator[](int i) const {
        return data[i];
    }

    letter& operator[](int i) {
        return data[i];
    }

    int size() const {
        return k;
    }

};

template<typename letter>
inline int compare(const kmer<letter> &s1, const kmer<letter> &s2) {
    for (int i = 0; ; i++) {
        if (i >= std::min(s1.size(), s2.size())) {
            return s1.size() - s2.size();
        }
        if (s1[i] != s2[i]) {
            return s1[i] - s2[i];
        }
    }
}



template<typename letter>
bool operator<(const kmer<letter> &s1, const kmer<letter> &s2) {
    return compare(s1, s2) < 0;
}

template<typename letter>
bool operator>(const kmer<letter> &s1, const kmer<letter> &s2) {
    return compare(s1, s2) > 0;
}

template<typename letter>
bool operator==(const kmer<letter> &s1, const kmer<letter> &s2) {
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

template<typename letter>
bool operator>=(const kmer<letter> &s1, const kmer<letter> &s2) {
    return !(s1 < s2);
}

template<typename letter>
bool operator<=(const kmer<letter> &s1, const kmer<letter> &s2) {
    return !(s1 > s2);
}

template<typename letter>
bool operator!=(const kmer<letter> &s1, const kmer<letter> &s2) {
    return !(s1 == s2);
}

namespace std {

    template <>
    template <typename letter>
    struct hash< kmer<letter> >
    {
        std::size_t operator()(const kmer<letter>& s) const {
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

template<typename letter>
struct kmer_count {
    kmer<letter> data;
    int count;

    kmer_count(const kmer<letter> &s, int c)
        : data(s)
        , count(c) {
    }

    kmer_count(const letter *seq = nullptr, int size = 0, int c = 0)
        : data(seq, size)
        , count(c) {
    }

    ~kmer_count() {
    }

};

template<typename letter>
vector1D< kmer_count<letter> > count_kmers(const vector1D<letter> &seq, int k) {

    tree_map<kmer<letter>, int> map;

    for (int i = 0; i < seq.size() - k + 1; i++) {
        kmer<letter> s(seq.data() + i, k);
        map[s]++;
    }

    int id = 0;
    vector1D< kmer_count<letter> > result(map.size());
    for (auto it = map.cbegin(); it != map.cend(); ++it) {
        result[id++] = kmer_count<letter>(it->first, it->second);
    }

    return result;
}

template<typename letter>
vector2D< kmer_count<letter> > count_all_kmers(const vector2D<letter> &seqs, int k) {

    vector2D< kmer_count<letter> > result(seqs.size(), 0);

    for (int i = 0; i < seqs.size(); i++) {
        result[i] = count_kmers(seqs[i], k);
    }

    return result;
}
