#include "format.hpp"

#include <set>
#include <map>

template<typename dtype>
sequence_array to_sequence_array(dtype *data, int rows, int cols) {

    std::set<int> symbols;
    for (int i = 0; i < rows * cols; i++) {
        symbols.insert(data[i]);
    }

    int symbol_size = 0;
    std::map<int, int> symbol_map;
    for (std::set<int>::iterator it = symbols.begin(); it != symbols.end(); ++it) {
        symbol_map[*it] = symbol_size++;
    }

    sequence_array seq_array = {
        symbol_size,
        cols,
        vector2D<ltype>(rows, cols)
    };

    vector2D<ltype> &sequences = seq_array.sequences;

    sequences.resize(rows);
    for (int i = 0; i < rows; i++) {
        sequences[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            sequences[i][j] = symbol_map[data[i * cols + j]];
        }
    }

    return seq_array;

}

template sequence_array to_sequence_array(int *, int, int);
