#include "format.hpp"

#include <set>
#include <map>

template<typename letter, typename dtype>
sequence_array<letter> to_sequence_array(dtype *data, int rows, int cols) {

    std::set<int> symbols;
    for (int i = 0; i < rows * cols; i++) {
        symbols.insert(data[i]);
    }

    int symbol_size = 0;
    std::map<int, int> symbol_map;
    for (std::set<int>::iterator it = symbols.begin(); it != symbols.end(); ++it) {
        symbol_map[*it] = symbol_size++;
    }

    sequence_array<letter> seq_array = {
        symbol_size,
        cols,
        vector2D<letter>(rows, cols)
    };
    // seq_array.alphabet_size =
    // seq_array.sequences_len = cols;

    vector2D<letter> &sequences = seq_array.sequences;

    sequences.resize(rows);
    for (int i = 0; i < rows; i++) {
        sequences[i].resize(cols);
        for (int j = 0; j < cols; j++) {
            sequences[i][j] = symbol_map[data[i * cols + j]];
        }
    }

    // sequence_array<letter> seq_array = {
    //     symbol_size,
    //     sequences[0].size(),
    //     std::move(sequences)
    // };

    return seq_array;

}

template sequence_array<int> to_sequence_array(int *, int, int);
