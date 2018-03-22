#include <cstring>
#include <algorithm>

#include "local_alignment.hpp"

struct dynamic_cell {
    float M;
    float X;
    float Y;
    // float X2;
    // float Y2;
};

template<typename T>
const T& max4(const T&a, const T&b, const T&c, const T&d) {
    using std::max;
    return max(max(a, b), max(c, d));
}

template<typename T>
const T& max3(const T&a, const T&b, const T&c) {
    using std::max;
    return max(a, max(b, c));
}

template<typename T>
const T& max2(const T&a, const T&b) {
    using std::max;
    return max(a, b);
}

template<typename dtype>
inline dtype local_alignment_seq(const vector1D<ltype> &sequence1, const vector1D<ltype> &sequence2,
                                 int alphabet_size, dtype substitution, dtype gap_open, dtype gap_extension) {

    matrix<dynamic_cell> table(sequence1.size() + 1, sequence2.size() + 1);
    // dynamic_cell *table = new dynamic_cell[(sequences_len + 1) * (sequences_len + 1)];

    for (int j = 0; j < sequence2.size() + 1; j++) {
        table(0, i) = {0.0f, 0.0f, 0.0f};
    }

    for (int i = 0; i < sequence1.size() + 1; i++) {
        table(i, 0) = {0.0f, 0.0f, 0.0f};
    }

    const dtype zero = dtype();

    for (int i = 0; i < sequence1.size(); i++) {

        for (int j = 0; j < sequence2.size(); j++) {


            dynamic_cell &cell = table(i + 1, j + 1);

            dtype cost = (sequence1[i] == sequence2[j]) ? substitution : -substitution;

            cell.M = max4(table(i, j).M + cost,
                          table(i, j).X + cost,
                          table(i, j).Y + cost,
                          zero);

            cell.X = max3(table(i, j + 1).X - gap_extension,
                          table(i, j + 1).M - gap_open,
                          zero);

            cell.Y = max3(table(i + 1, j).Y - gap_extension,
                          table(i + 1, j).M - gap_open,
                          zero);

        }

    }

    const dynamic_cell &result = table(sequence1.size(), sequence2.size());

    return max3(result.M, result.X, result.Y);

}



template<typename dtype>
sq_matrix<dtype> local_alignment(const vector2D<ltype> &sequences,
                                 int sequences_len, int alphabet_size,
                                 dtype substitution, dtype gap_open,
                                 dtype gap_extension) {

    sq_matrix<dtype> K(sequences.size(), sequences.size());

    for (int i = 0; i < sequences.size(); i++) {

        for (int j = i; j < sequences.size(); j++) {

            K(i, j) = local_alignment_seq(sequences[i], sequences[j], alphabet_size,
                                          substitution, gap_open, gap_extension);

        }
    }

    K.symmetrise_lower();

    return K;
}


template sq_matrix<double> local_alignment(const vector2D<ltype> &, int, int, double, double, double);
template sq_matrix<float> local_alignment(const vector2D<ltype> &, int, int, float, float, float);
