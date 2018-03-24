#include <cstring>
#include <algorithm>
#include <limits>

#include "local_alignment.hpp"


template<typename T>
const T& max(const T&a, const T&b, const T&c, const T&d) {
    using std::max;
    return max(max(a, b), max(c, d));
}

template<typename T>
const T& max(const T&a, const T&b, const T&c) {
    using std::max;
    return max(a, max(b, c));
}

template<typename T>
const T& max(const T&a, const T&b) {
    using std::max;
    return max(a, b);
}

template<typename dtype>
struct dynamic_cell {
    dtype M;
    dtype X;
    dtype Y;
};

template<typename dtype>
inline dtype needleman_wunsch_affine_gap_pair(const vector1D<ltype> &sequence1,
                                              const vector1D<ltype> &sequence2,
                                              int alphabet_size,
                                              dtype similarity, dtype substitution,
                                              dtype gap_open, dtype gap_extension) {

    using cell_type = dynamic_cell<float>;

    // TODO(RL) Do something !
    // if (sequence1.size() > sequence2.size()) {
    //     std::swap(sequence1, sequence2);
    // }

    vector1D<cell_type> table1(sequence2.size() + 1, sequence2.size() + 1);
    vector1D<cell_type> table2(sequence2.size() + 1, sequence2.size() + 1);

    const dtype zero = dtype();
    const dtype infm = std::numeric_limits<dtype>::min();
    // const dtype infp = std::numeric_limits<dtype>::max();

    table1[0] = {zero, zero, zero};
    for (int j = 0; j < sequence2.size(); j++) {
        table1[j + 1] = {infm, - (j + 1) * gap_extension, infm};
    }

    for (int i = 0; i < sequence1.size(); i++) {

        table2[0] = {infm, infm, - (i + 1) * gap_extension};

        for (int j = 0; j < sequence2.size(); j++) {

            cell_type &cell = table2[j + 1];

            dtype cost = (sequence1[i] == sequence2[j]) ? similarity : -substitution;

            // TODO(RL) This formula I've derived myself (verify)
            cell.M = max(table1[j].M + cost,
                         table1[j].X + cost,
                         table1[j].Y + cost);

            cell.X = max(table1[j+1].X - gap_extension,
                         table1[j+1].M - gap_open);

            cell.Y = max(table2[j].Y - gap_extension,
                         table2[j].M - gap_open);


        }

        std::swap(table1, table2);

    }

    cell_type &result = table1[sequence2.size()];

    return max(result.M, result.X, result.Y);

}

// TODO(RL) Pass a similarity matrix directly
template<typename dtype>
inline dtype smith_waterman_affine_gap_pair(const vector1D<ltype> &sequence1,
                                            const vector1D<ltype> &sequence2,
                                            int alphabet_size,
                                            dtype similarity, dtype substitution,
                                            dtype gap_open, dtype gap_extension) {

    using cell_type = dynamic_cell<float>;

    // TODO(RL) Do something !
    // if (sequence1.size() > sequence2.size()) {
    //     std::swap(sequence1, sequence2);
    // }

    vector1D<cell_type> table1(sequence2.size() + 1, sequence2.size() + 1);
    vector1D<cell_type> table2(sequence2.size() + 1, sequence2.size() + 1);

    const dtype zero = dtype();

    for (int j = 0; j < sequence2.size() + 1; j++) {
        table1[j] = {zero, zero, zero};
    }

    table2[0] = {zero, zero, zero};

    dtype best_score = zero;
    for (int i = 0; i < sequence1.size(); i++) {

        for (int j = 0; j < sequence2.size(); j++) {

            cell_type &cell = table2[j + 1];

            dtype cost = (sequence1[i] == sequence2[j]) ? similarity : -substitution;

            // TODO(RL) This formula I've derived myself (verify)
            cell.M = max(table1[j].M + cost,
                         table1[j].X + cost,
                         table1[j].Y + cost,
                         zero);

            cell.X = max(table1[j+1].X - gap_extension,
                         table1[j+1].M - gap_open);

            cell.Y = max(table2[j].Y - gap_extension,
                         table2[j].M - gap_open);

            best_score = max(best_score, cell.M, cell.X, cell.Y);

        }

        std::swap(table1, table2);

    }

    return best_score;

}

template<typename dtype>
sq_matrix<dtype> smith_waterman_affine_gap(const vector2D<ltype> &sequences,
                                           int sequences_len, int alphabet_size,
                                           dtype similarity, dtype substitution,
                                           dtype gap_open, dtype gap_extension) {

    sq_matrix<dtype> K(sequences.size(), sequences.size());

    for (int i = 0; i < sequences.size(); i++) {

        for (int j = i; j < sequences.size(); j++) {

            K(i, j) = smith_waterman_affine_gap_pair(sequences[i], sequences[j],
                                                     alphabet_size,
                                                     similarity, substitution,
                                                     gap_open, gap_extension);
        }
    }

    K.symmetrise_lower();

    return K;
}


template<typename dtype>
sq_matrix<dtype> needleman_wunsch_affine_gap(const vector2D<ltype> &sequences,
                                             int sequences_len, int alphabet_size,
                                             dtype similarity, dtype substitution,
                                             dtype gap_open, dtype gap_extension) {

    sq_matrix<dtype> K(sequences.size(), sequences.size());

    for (int i = 0; i < sequences.size(); i++) {

        for (int j = i; j < sequences.size(); j++) {

            K(i, j) = needleman_wunsch_affine_gap_pair(sequences[i], sequences[j],
                                                       alphabet_size,
                                                       similarity, substitution,
                                                       gap_open, gap_extension);
        }
    }

    K.symmetrise_lower();

    return K;
}

// template sq_matrix<double> local_alignment(const vector2D<ltype> &, int, int, double, double, double);
template sq_matrix<float> needleman_wunsch_affine_gap(const vector2D<ltype> &, int, int, float, float, float, float);
template sq_matrix<float> smith_waterman_affine_gap(const vector2D<ltype> &, int, int, float, float, float, float);
