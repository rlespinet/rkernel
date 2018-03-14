#include "spectrum.hpp"
#include "utils.hpp"

#include "kmer.hpp"

template<typename letter, typename dtype>
sq_matrix<dtype> spectrum(const vector2D<letter> &sequences,
                   int sequences_len, int alphabet_size, int k) {

    vector2D< kmer_count<letter> > kmers = count_all_kmers(sequences, k);

    sq_matrix<dtype> K(sequences.size());

    // Note(RL) : Below is the critical part

    for (int i = 0; i < kmers.size(); i++) {

        for (int j = i; j < kmers.size(); j++) {

            int matches = 0;
            int ki = 0, kj = 0;
            while (ki < kmers[i].size() && kj < kmers[j].size()) {
                const int cmp = compare(kmers[i][ki].data, kmers[j][kj].data);
                if (cmp < 0) {
                    ki++;
                } else if (cmp > 0) {
                    kj++;
                } else {
                    matches += kmers[i][ki].count * kmers[j][kj].count;
                    ki++;
                    kj++;
                }
            }

            K(i, j) = matches;

        }

    }

    K.symmetrise_lower();

    return K;
}


// template<typename letter, typename dtype>
// matrix<dtype> spectrum(const vector2D<letter> &sequences1,
//                        const vector2D<letter> &sequences2,
//                        int sequences_len, int alphabet_size, int k) {

//     vector2D<kmer<letter> > kmers1 = kmer_encode_all(sequences1, k, sequences_len, alphabet_size);
//     vector2D<kmer<letter> > kmers2 = kmer_encode_all(sequences2, k, sequences_len, alphabet_size);

//     matrix<dtype> K(sequences1.size(), sequences2.size());

//     // Note(RL) : Below is the critical part

//     for (int i = 0; i < kmers1.size(); i++) {

//         for (int j = 0; j < kmers2.size(); j++) {

//             int matches = 0;
//             int ki = 0, kj = 0;
//             while (ki < kmers1[i].size() && kj < kmers2[j].size()) {
//                 const int cmp = compare(kmers1[i][ki], kmers2[j][kj]);
//                 if (cmp < 0) {
//                     ki++;
//                 } else if (cmp > 0) {
//                     kj++;
//                 } else {
//                     matches += kmers1[i][ki].count * kmers2[j][kj].count;
//                     ki++;
//                     kj++;
//                 }
//             }

//             K(i, j) = matches;

//         }

//     }

//     return K;
// }

// Template specialization

template sq_matrix<double> spectrum(const vector2D<int> &, int, int, int);
template sq_matrix<float> spectrum(const vector2D<int> &, int, int, int);

// template matrix<double> spectrum(const vector2D<int> &, const vector2D<int> &, int, int, int);
// template matrix<float> spectrum(const vector2D<int> &, const vector2D<int> &, int, int, int);
