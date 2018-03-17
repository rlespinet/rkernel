#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "wildcard.hpp"

struct kmer_wildcard : kmer_count {
    int seq_id;

    kmer_wildcard()
	: kmer_count()
	, seq_id(0) {
    }

    kmer_wildcard(const kmer_wildcard &oth)
        : kmer_count(oth)
        , seq_id(oth.seq_id) {
    }

    kmer_wildcard(const kmer_count &s, int i)
	: kmer_count(s)
	, seq_id(i) {
    }

};

inline bool compare(const kmer_wildcard &s1, const kmer_wildcard &s2) {
    return s1.seq_id < s2.seq_id;
}

vector1D< kmer_wildcard > compute_kmer_wildcard(const vector2D<ltype> &sequences,
                                                int sequences_len, int alphabet_size,
                                                int k) {

    vector2D< kmer_count > kmers = count_all_kmers(sequences, k);

    int total_kmers = 0;
    for (int i = 0; i < kmers.size(); i++) {
        total_kmers += kmers[i].size();
    }

    vector1D< kmer_wildcard > kmer_wildcards(total_kmers);
    for (int i = 0; i < kmers.size(); i++) {
        for (int j = 0; j < kmers[i].size(); j++) {
            const kmer_wildcard elt(kmers[i][j], i);
            kmer_wildcards.push_back(elt);
        }
    }

    return kmer_wildcards;
}

template<typename dtype>
void wildcard_compute_rec(sq_matrix<dtype> &K,
                          const vector1D< kmer_wildcard > &tracks,
                          int alphabet_size, int k, int d, int m,
                          float w, float wc = 1.0) {

    if (tracks.size() == 0) {
        return;
    }

    if (d == k) {
        // NOTE(RL) Normally, we would iterate for i=0..size and j=0..size
        // but in fact we can iterator for j=i..size, and multiply
        // by a coefficient (coeff) to compensate
        for (int i = 0; i < tracks.size(); i++) {
            for (int j = i; j < tracks.size(); j++) {

                int id1 = tracks[i].seq_id;
                int id2 = tracks[j].seq_id;

                int matches = tracks[i].count * tracks[j].count;

                float coeff = (i != j && id1 == id2) ? 2.0f : 1.0f;
                K(id1, id2) += coeff * matches * wc * wc;
            }
        }

        return;
    }

    if (m > 0) {
        wildcard_compute_rec<dtype>(K, tracks, alphabet_size, k, d+1, m-1, w, wc * w);
    }

    vector1D< kmer_wildcard > new_tracks(tracks.size());

    for (ltype a = 0; a < alphabet_size; a++) {

        new_tracks.clear();

        for (int i = 0; i < tracks.size(); i++) {

            ltype c = tracks[i].data[d];
            if (c == a) {
                new_tracks.push_back(tracks[i]);
            }

        }

        wildcard_compute_rec<dtype>(K, new_tracks, alphabet_size, k, d+1, m, w, wc);

    }

}


template<typename dtype>
sq_matrix<dtype> wildcard(const vector2D<ltype> &sequences,
                          int sequences_len, int alphabet_size,
                          int k, int m, float w) {

    vector1D<kmer_wildcard > kmer_wildcards = compute_kmer_wildcard(sequences, sequences_len,
                                                                    alphabet_size, k);
    sq_matrix<dtype> K(sequences.size(), 0);

    wildcard_compute_rec<dtype>(K, kmer_wildcards, alphabet_size, k, 0, m, w);

    K.symmetrise_lower();

    return K;

}


template sq_matrix<double> wildcard(const vector2D<ltype> &, int, int, int, int, float);
template sq_matrix<float> wildcard(const vector2D<ltype> &, int, int, int, int, float);
