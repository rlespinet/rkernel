#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "mismatch.hpp"

struct kmer_mismatch : kmer {
    using kmer::encoding;
    using kmer::count;
    int seq_id;
    int mismatchs;

    kmer_mismatch(kmer k = kmer(), int s = 0, int m = 0)
        : kmer(k)
        , seq_id(s)
        , mismatchs(m) {}

};

inline bool compare(const kmer_mismatch &s1, const kmer_mismatch &s2) {
    return s1.seq_id < s2.seq_id;
}

template<typename letter>
vector1D<kmer_mismatch> compute_kmer_mismatch(const vector2D<letter> &sequences,
                                              int sequences_len, int alphabet_size,
                                              int k, int m) {

    vector2D<kmer> kmers = kmer_encode_all(sequences, k, sequences_len, alphabet_size);

    int total_kmers = 0;
    for (int i = 0; i < kmers.size(); i++) {
        total_kmers += kmers[i].size();
    }

    int kmer_mismatch_size = 0;
    vector1D<kmer_mismatch> kmer_mismatchs(total_kmers);
    for (int i = 0; i < kmers.size(); i++) {
        for (int j = 0; j < kmers[i].size(); j++) {
            kmer_mismatchs[kmer_mismatch_size] = kmer_mismatch(kmers[i][j], i, m);
            kmer_mismatch_size++;
        }
    }
    kmer_mismatchs.resize(kmer_mismatch_size);

    return kmer_mismatchs;
}

template<typename letter, typename dtype>
void mismatch_compute_rec(sq_matrix<dtype> &K,
                          const vector1D<kmer_mismatch> &tracks,
                          int alphabet_size, int k, int d) {

    if (tracks.size() == 0) {
        return;
    }

    if (d == k) {
        // This is a leaf !
        // std::sort(tracks, tracks + len, compare<hash_t>);

        // NOTE(RL) Normally, we would iterate for i=0..size and j=0..size
        // but in fact we can iterator for j=i..size, and multiply
        // by a coefficient (coeff) to compensate
        for (int i = 0; i < tracks.size(); i++) {
            for (int j = i; j < tracks.size(); j++) {

                int id1 = tracks[i].seq_id;
                int id2 = tracks[j].seq_id;

                int matches = tracks[i].count * tracks[j].count;

                int coeff = (i != j && id1 == id2) ? 2 : 1;
                K(id1, id2) += coeff * matches;
            }
        }

        return;
    }

    vector1D<kmer_mismatch> new_tracks(tracks.size());

    for (letter a = 0; a < alphabet_size; a++) {

        int new_tracks_size = 0;
        for (int i = 0; i < tracks.size(); i++) {

            letter c = kmer_decode<letter>(tracks[i].encoding, alphabet_size, k - 1 - d);
            if (c == a || tracks[i].mismatchs > 0) {

                new_tracks[new_tracks_size] = tracks[i];
                if (c != a) {
                    new_tracks[new_tracks_size].mismatchs --;
                }

                new_tracks_size++;
            }

        }

        new_tracks.resize(new_tracks_size);

        mismatch_compute_rec<letter, dtype>(K, new_tracks, alphabet_size, k, d + 1);

    }

}


template<typename letter, typename dtype>
sq_matrix<dtype> mismatch(const vector2D<letter> &sequences,
                          int sequences_len, int alphabet_size,
                          int k, int m) {

    vector1D<kmer_mismatch> kmer_mismatchs = compute_kmer_mismatch(sequences, sequences_len,
                                                                   alphabet_size, k, m);
    sq_matrix<dtype> K(sequences.size(), 0);

    mismatch_compute_rec<letter, dtype>(K, kmer_mismatchs, alphabet_size, k, 0);

    K.symmetrise_lower();

    return K;

}


template sq_matrix<double> mismatch(const vector2D<int> &, int, int, int, int);
template sq_matrix<float> mismatch(const vector2D<int> &, int, int, int, int);

// template matrix<double> spectrum(const vector2D<int> &, const vector2D<int> &, int, int, int);
// template matrix<float> spectrum(const vector2D<int> &, const vector2D<int> &, int, int, int);
