#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "dimismatch.hpp"

struct kmer_dimismatch : kmer_count {
    int seq_id;
    int mismatchs;
    bool last;
    int total;

    kmer_dimismatch()
	: kmer_count()
	, seq_id(0)
	, mismatchs(0)
        , last(false)
        , total(0) {
    }

    kmer_dimismatch(const kmer_dimismatch &oth)
        : kmer_count(oth)
        , seq_id(oth.seq_id)
        , mismatchs(oth.mismatchs)
        , last(oth.last)
        , total(oth.total) {
    }

    kmer_dimismatch(const kmer_count &s, int i, int m, bool last = false, int total = 0)
	: kmer_count(s)
	, seq_id(i)
	, mismatchs(m)
        , last(last)
        , total(total) {
    }

};

inline bool compare(const kmer_dimismatch &s1, const kmer_dimismatch &s2) {
    return s1.seq_id < s2.seq_id;
}

vector1D< kmer_dimismatch > compute_kmer_dimismatch(const vector2D<ltype> &sequences,
                                                  int sequences_len, int alphabet_size,
                                                  int k, int m) {

    vector2D< kmer_count > kmers = count_all_kmers(sequences, k);

    int total_kmers = 0;
    for (int i = 0; i < kmers.size(); i++) {
        total_kmers += kmers[i].size();
    }

    vector1D< kmer_dimismatch > kmer_dimismatchs(total_kmers);
    for (int i = 0; i < kmers.size(); i++) {
        for (int j = 0; j < kmers[i].size(); j++) {
            const kmer_dimismatch elt(kmers[i][j], i, m);
            kmer_dimismatchs.push_back(elt);
        }
    }

    return kmer_dimismatchs;
}


template<typename dtype>
void dimismatch_compute_rec(sq_matrix<dtype> &K,
                          const vector1D< kmer_dimismatch > &tracks,
                          int alphabet_size, int k, int d) {

    if (tracks.size() == 0) {
        return;
    }

    if (d == k) {

        vector1D<int> matches(tracks.size());
        vector1D<int> ids(tracks.size());

        matches.push_back(0);
        ids.push_back(tracks[0].seq_id);
        for (int i = 0; i < tracks.size(); i++) {

            int id = tracks[i].seq_id;

            if (id != ids.last()) {
                matches.push_back(0);
                ids.push_back(id);
            }
            matches.last() += tracks[i].count * tracks[i].total;
        }

        for (int i = 0; i < matches.size(); i++) {

            for (int j = i; j < matches.size(); j++) {

                int idi = ids[i];
                int idj = ids[j];

                K(idi, idj) += matches[i] * matches[j];
            }
        }

        return;
    }

    vector1D< kmer_dimismatch > new_tracks(tracks.size());

    for (ltype a = 0; a < alphabet_size; a++) {

        new_tracks.clear();

        for (int i = 0; i < tracks.size(); i++) {

            ltype c = tracks[i].data[d];
            if (c == a || tracks[i].mismatchs > 0) {

                int new_mismatchs = tracks[i].mismatchs - ((c == a) ? 0 : 1);
                bool new_last = (c == a);
                int new_total = tracks[i].total + (tracks[i].last && new_last ? 1 : 0);

                const kmer_dimismatch new_kmer(tracks[i], tracks[i].seq_id,
                                               new_mismatchs, new_last, new_total);

                new_tracks.push_back(new_kmer);

            }

        }

        dimismatch_compute_rec<dtype>(K, new_tracks, alphabet_size, k, d + 1);

    }

}

template<typename dtype>
sq_matrix<dtype> dimismatch(const vector2D<ltype> &sequences,
                            int sequences_len, int alphabet_size,
                            int k, int m) {

    vector1D<kmer_dimismatch > kmer_dimismatchs = compute_kmer_dimismatch(sequences, sequences_len,
                                                                          alphabet_size, k, m);
    sq_matrix<dtype> K(sequences.size(), 0);

    dimismatch_compute_rec<dtype>(K, kmer_dimismatchs, alphabet_size, k, 0);

    K.symmetrise_lower();

    return K;

}

template sq_matrix<double> dimismatch(const vector2D<ltype> &, int, int, int, int);
template sq_matrix<float> dimismatch(const vector2D<ltype> &, int, int, int, int);
