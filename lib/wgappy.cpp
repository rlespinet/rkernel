#include <cstring>
#include <algorithm>

#include "kmer.hpp"
#include "wgappy.hpp"

struct kmer_wgappy : kmer_count {
    int seq_id;
    int pos_id;

    kmer_wgappy()
	: kmer_count()
	, seq_id(0)
	, pos_id(0) {
    }

    kmer_wgappy(const kmer_count &s, int i, int p = 0)
	: kmer_count(s)
	, seq_id(i)
	, pos_id(p) {
    }

};

inline bool compare_wgappy(const kmer_wgappy &s1, const kmer_wgappy &s2) {
    return s1.seq_id < s2.seq_id;
}

inline float substring_occurences(const kmer &str, const kmer &substr, float weight = 1.0f) {

    const int buffer_size = str.size() - substr.size() + 1;
    vector1D<float> buffer(buffer_size, buffer_size);

    for (int j = 0; j < buffer.size(); j++) {
        buffer[j] = 1.0f;
    }

    for (int i = 0; i < substr.size(); i++) {

        buffer[0] = (str[i] == substr[i]) ? buffer[0] : 0.0f;
        const float coeff = (i < substr.size() - 1) ? weight : 1.0f;

        for (int j = 1; j < buffer.size(); j++) {
            const float match = (str[i + j] == substr[i]) ? buffer[j] : 0.0f;
            buffer[j] = buffer[j-1] * coeff + match;
        }
    }

    return buffer[buffer_size - 1];

}

vector1D< kmer_wgappy > compute_kmer_wgappy(const vector2D<ltype> &sequences,
                                            int sequences_len, int alphabet_size,
                                            int g) {

    vector2D< kmer_count > kmers = count_all_kmers(sequences, g);

    int total_kmers = 0;
    for (int i = 0; i < kmers.size(); i++) {
        total_kmers += kmers[i].size();
    }

    vector1D< kmer_wgappy > kmer_wgappys(total_kmers);
    for (int i = 0; i < kmers.size(); i++) {
        for (int j = 0; j < kmers[i].size(); j++) {
            const kmer_wgappy elt(kmers[i][j], i);
            kmer_wgappys.push_back(elt);
        }
    }

    return kmer_wgappys;
}

template<typename dtype>
void wgappy_compute_rec(sq_matrix<dtype> &K,
                        vector1D< kmer_wgappy > &tracks,
                        vector1D<ltype> &branch,
                        int alphabet_size, int k, float w) {

    if (tracks.size() == 0) {
        return;
    }

    if (branch.size() == k) {

        kmer cmp(branch.data(), branch.size());

        vector1D<float> matches(tracks.size());
        vector1D<int> ids(tracks.size());

        matches.push_back(0.0f);
        ids.push_back(tracks[0].seq_id);
        for (int i = 0; i < tracks.size(); i++) {

            int id = tracks[i].seq_id;

            if (id != ids.last()) {
                matches.push_back(0.0f);
                ids.push_back(id);
            }
            matches.last() += substring_occurences(tracks[i].data, cmp, w) * tracks[i].count;
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

    vector1D< kmer_wgappy > new_tracks(tracks.size());

    for (ltype a = 0; a < alphabet_size; a++) {

        new_tracks.clear();

        for (int i = 0; i < tracks.size(); i++) {

            int next_a = find(tracks[i].data, a, tracks[i].pos_id);
            if (next_a < 0) {
                continue;
            }

            const kmer_wgappy track(tracks[i], tracks[i].seq_id, next_a + 1);
            new_tracks.push_back(track);
        }

        branch.push_back(a);
        wgappy_compute_rec<dtype>(K, new_tracks, branch, alphabet_size, k, w);
        branch.pop_back();

    }

}


template<typename dtype>
sq_matrix<dtype> wgappy(const vector2D<ltype> &sequences,
                        int sequences_len, int alphabet_size,
                        int k, int g, float w) {

    vector1D<kmer_wgappy > kmer_wgappys = compute_kmer_wgappy(sequences, sequences_len,
                                                              alphabet_size, g);
    sq_matrix<dtype> K(sequences.size(), 0);

    vector1D<ltype> branch(k);
    wgappy_compute_rec<dtype>(K, kmer_wgappys, branch, alphabet_size, k, w);

    K.symmetrise_lower();

    return K;

}


template sq_matrix<double> wgappy(const vector2D<ltype> &, int, int, int, int, float);
template sq_matrix<float> wgappy(const vector2D<ltype> &, int, int, int, int, float);
