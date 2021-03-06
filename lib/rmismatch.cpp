#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "rmismatch.hpp"

struct kmer_rmismatch : kmer_count {
    int seq_id;
    int mismatchs;
    int consecutive;
    int total;

    kmer_rmismatch()
	: kmer_count()
	, seq_id(0)
	, mismatchs(0)
        , consecutive(0)
        , total(0) {
    }

    kmer_rmismatch(const kmer_rmismatch &oth)
        : kmer_count(oth)
        , seq_id(oth.seq_id)
        , mismatchs(oth.mismatchs)
        , consecutive(oth.consecutive)
        , total(oth.total) {
    }

    kmer_rmismatch(const kmer_count &s, int i, int m, int consecutive = 0, int total = 0)
	: kmer_count(s)
	, seq_id(i)
	, mismatchs(m)
        , consecutive(consecutive)
        , total(total) {
    }

};

inline bool compare(const kmer_rmismatch &s1, const kmer_rmismatch &s2) {
    return s1.seq_id < s2.seq_id;
}

vector1D< kmer_rmismatch > compute_kmer_rmismatch(const vector2D<ltype> &sequences,
                                                  int sequences_len, int alphabet_size,
                                                  int k, int m) {

    vector2D< kmer_count > kmers = count_all_kmers(sequences, k);

    int total_kmers = 0;
    for (int i = 0; i < kmers.size(); i++) {
        total_kmers += kmers[i].size();
    }

    vector1D< kmer_rmismatch > kmer_rmismatchs(total_kmers);
    for (int i = 0; i < kmers.size(); i++) {
        for (int j = 0; j < kmers[i].size(); j++) {
            const kmer_rmismatch elt(kmers[i][j], i, m);
            kmer_rmismatchs.push_back(elt);
        }
    }

    return kmer_rmismatchs;
}


template<typename dtype>
void rmismatch_compute_rec(sq_matrix<dtype> &K,
                           const vector1D< kmer_rmismatch > &tracks,
                           int alphabet_size, int k, int r, int d) {

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

    vector1D< kmer_rmismatch > new_tracks(tracks.size());

    for (ltype a = 0; a < alphabet_size; a++) {

        new_tracks.clear();

        for (int i = 0; i < tracks.size(); i++) {

            ltype c = tracks[i].data[d];
            if (c == a || tracks[i].mismatchs > 0) {

                int new_mismatchs = tracks[i].mismatchs - ((c == a) ? 0 : 1);
                int new_consecutive = (c == a) ? tracks[i].consecutive + 1 : 0;
                int new_total = tracks[i].total + (new_consecutive >= r ? 1 : 0);

                const kmer_rmismatch new_kmer(tracks[i], tracks[i].seq_id,
                                              new_mismatchs, new_consecutive, new_total);

                new_tracks.push_back(new_kmer);

            }

        }

        rmismatch_compute_rec<dtype>(K, new_tracks, alphabet_size, k, r, d + 1);

    }

}

template<typename dtype>
sq_matrix<dtype> rmismatch(const vector2D<ltype> &sequences,
                           int sequences_len, int alphabet_size,
                           int k, int m, int r) {

    vector1D<kmer_rmismatch > kmer_rmismatchs = compute_kmer_rmismatch(sequences, sequences_len,
                                                                       alphabet_size, k, m);
    sq_matrix<dtype> K(sequences.size(), 0);

    rmismatch_compute_rec<dtype>(K, kmer_rmismatchs, alphabet_size, k, r, 0);

    K.symmetrise_lower();

    return K;

}

template<typename dtype>
sq_matrix<dtype> dimismatch(const vector2D<ltype> &sequences,
                            int sequences_len, int alphabet_size,
                            int k, int m) {
    return rmismatch<dtype>(sequences, sequences_len, alphabet_size, k, m, 2);
}

template sq_matrix<double> rmismatch(const vector2D<ltype> &, int, int, int, int, int);
template sq_matrix<float> rmismatch(const vector2D<ltype> &, int, int, int, int, int);

template sq_matrix<double> dimismatch(const vector2D<ltype> &, int, int, int, int);
template sq_matrix<float> dimismatch(const vector2D<ltype> &, int, int, int, int);
