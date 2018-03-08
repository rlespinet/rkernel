#include "format.hpp"

#include <set>
#include <map>

seq_array_t format_sequence_array(int *data, int rows, int cols) {

    std::set<int> symbols;
    for (int i = 0; i < rows * cols; i++) {
        symbols.insert(data[i]);
    }

    int count = 0;
    std::map<int, int> symbol_map;
    for (std::set<int>::iterator it = symbols.begin(); it != symbols.end(); ++it) {
        symbol_map[*it] = count++;
    }

    int *reindexed_data = new int[rows * cols];
    if (reindexed_data == NULL) {
        return invalid_seq_array;
    }

    for (int i = 0; i < rows * cols; i++) {
        reindexed_data[i] = symbol_map[data[i]];
    }

    seq_array_t seq_array = {
        reindexed_data,
        rows, cols,
        count
    };

    return seq_array;

}
