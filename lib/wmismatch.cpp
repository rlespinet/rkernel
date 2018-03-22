#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "wmismatch.hpp"

struct kmer_mismatch : kmer_count {
    int seq_id;
    int mismatchs;

    kmer_mismatch()
	: kmer_count()
	, seq_id(0)
	, mismatchs(0) {
    }

    kmer_mismatch(const kmer_mismatch &oth)
        : kmer_count(oth)
        , seq_id(oth.seq_id)
        , mismatchs(oth.mismatchs) {
    }

    kmer_mismatch(const kmer_count &s, int i, int m)
	: kmer_count(s)
	, seq_id(i)
	, mismatchs(m) {
    }

};

inline bool compare(const kmer_mismatch &s1, const kmer_mismatch &s2) {
    return s1.seq_id < s2.seq_id;
}

vector1D< kmer_mismatch > compute_kmer_mismatch(const vector2D<ltype> &sequences,
                                                int sequences_len, int alphabet_size,
                                                int k) {

    vector2D< kmer_count > kmers = count_all_kmers(sequences, k);

    int total_kmers = 0;
    for (int i = 0; i < kmers.size(); i++) {
        total_kmers += kmers[i].size();
    }

    vector1D< kmer_mismatch > kmer_mismatchs(total_kmers);
    for (int i = 0; i < kmers.size(); i++) {
        for (int j = 0; j < kmers[i].size(); j++) {
            const kmer_mismatch elt(kmers[i][j], i, 0);
            kmer_mismatchs.push_back(elt);
        }
    }

    return kmer_mismatchs;
}

template<typename dtype>
void wmismatch_compute_rec(sq_matrix<dtype> &K,
                           const vector1D< kmer_mismatch > &tracks,
                           int alphabet_size, int k,
                           const float* weights, int m, int d) {

    if (tracks.size() == 0) {
        return;
    }

    if (d == k) {

        int *matches = new int[(m + 1) * tracks.size()];
        int *ids = new int[tracks.size()];

        int off_ids = 0;
        int off_matches = 0;

        std::fill(matches, matches + (m+1), 0);
        ids[0] = tracks[0].seq_id;

        for (int i = 0; i < tracks.size(); i++) {

            int id = tracks[i].seq_id;

            if (id != ids[off_ids]) {
                off_ids += 1;
                off_matches += m + 1;
                ids[off_ids] = id;
                std::fill(matches + off_matches, matches + off_matches + m + 1, 0);
            }

            matches[off_matches + tracks[i].mismatchs] += tracks[i].count;
        }

        for (int i = 0; i < off_ids + 1; i++) {

            for (int j = i; j < off_ids + 1; j++) {

                int idi = ids[i];
                int idj = ids[j];

                const int *matches_i = matches + (m+1) * i;
                const int *matches_j = matches + (m+1) * j;


                for (int k = 0; k < m + 1; k++) {
                    K(idi, idj) += (weights[k] * matches_i[k] *
                                    weights[k] * matches_j[k]);

                }

            }
        }

        delete[] matches;
        delete[] ids;

        return;
    }

    vector1D< kmer_mismatch > new_tracks(tracks.size());

    for (ltype a = 0; a < alphabet_size; a++) {

        new_tracks.clear();

        for (int i = 0; i < tracks.size(); i++) {

            ltype c = tracks[i].data[d];
            if (c == a || tracks[i].mismatchs < m) {

                new_tracks.push_back(tracks[i]);
                if (c != a) {
                    new_tracks.last().mismatchs++;
                }

            }

        }

        wmismatch_compute_rec<dtype>(K, new_tracks, alphabet_size, k, weights, m, d + 1);

    }

}


template<typename dtype>
sq_matrix<dtype> wmismatch(const vector2D<ltype> &sequences,
                           int sequences_len, int alphabet_size,
                           int k, int m, float weight) {

    vector1D<kmer_mismatch > kmer_mismatchs = compute_kmer_mismatch(sequences, sequences_len,
                                                                    alphabet_size, k);
    sq_matrix<dtype> K(sequences.size(), 0);

    float *weights = new float[m + 1];

    weights[0] = 1.0f;
    for (int k = 1; k < m + 1; k++) {
        weights[k] = weights[k-1] * weight;
    }

    wmismatch_compute_rec<dtype>(K, kmer_mismatchs, alphabet_size, k, weights, m, 0);

    K.symmetrise_lower();

    return K;

}


template sq_matrix<double> wmismatch(const vector2D<ltype> &, int, int, int, int, float);
template sq_matrix<float> wmismatch(const vector2D<ltype> &, int, int, int, int, float);
