#include <iostream>
#include <limits>
#include <vector>
#include <cassert>

#include "tree.hpp"

tree<int> construct_tree(std::vector<int> &s, int k, int l) {
    tree<int> t(l, k);

    for (int i = 0; i < s.size() - k + 1; i++) {

        tree<int>::node node = t.root();
        for (int j = i; j < i + k; j++) {
            node = t.child(node, s[j]);
        }
        t[node]++;

    }

    return t;
}

int compare_to_tree_rec(tree<int> &t, tree<int>::node node, const int *data, int len, int m) {

    if (m < 0) {
        return 0;
    }

    if (len == 0 || t.is_leaf(node)) {
        // In fact these 2 conditions should be true at the same time
        assert(len == 0 && t.is_leaf(node));
        return 0;
    }

    int matches = 0;

    for (int i = 0; i < t.N; i++) {

        tree<int>::node child = t.child(node, i);

        int m_child = (i == data[0]) ? m : m - 1;
        matches += compare_to_tree_rec(t, child, data + 1, len - 1, m_child) + t[child];

    }

    return matches;

}

int compare_to_tree(tree<int> &t, std::vector<int> &s, int m) {

    const int *data = s.data();

    int matches = 0;
    for (int i = 0; i < s.size() - t.D + 1; i++) {
        matches += compare_to_tree_rec(t, t.root(), data + i, t.D, m);
    }

    return matches;

}

std::vector<int> convert_dna(const char* s, const int *mapping) {

    std::vector<int> result;

    for (const char* p = s;*p;p++) {
        result.push_back(mapping[static_cast<int>(*p)]);
    }

    return result;

}

// void debug(tree<int> &t, tree<int>::node n, std::vector<int> &s) {

//     if (t.is_leaf(n)) {
//         return;
//     }

//     for (int i = 0; i < t.N; i++) {

//         s.push_back(i);
//         tree<int>::node m = t.child(n, i);
//         if (t[m] != 0) {
//             for (int j = 0; j < s.size(); j++) {
//                 std::cout << s[j] << " ";
//             }
//             std::cout << std::endl;
//         }

//         debug(t, m, s);
//         s.pop_back();
//     }
// }

int *compute_kernel(const char**strs, int len) {

    int *mapping = new int[256];

    mapping['A'] = 0;
    mapping['T'] = 1;
    mapping['G'] = 2;
    mapping['C'] = 3;

    strlen
    for (int i = 1; i < len; i++) {

    }


    // TODO(RL) allocate one bloc !
    int **x = new int*[len];
    for (int i = 0; i < len; i++) {
        x[i] = new int[len];
    }

    for (int i = 0; i < len; i++) {


    }


    for (int i = 0; i < len; i++) {
        delete[] x[i];
    }

    delete[] x;

    delete[] mapping;

}


int main(int argc, char**argv) {

    int *mapping = new int[256];

    // Initialize to INT_MIN, hopefully it will make
    // the program crash if something goes wrong
    for (int i = 0; i < 256; i++) {
        mapping[i] = std::numeric_limits<int>::min();
    }

    mapping['A'] = 0;
    mapping['T'] = 1;
    mapping['G'] = 2;
    mapping['C'] = 3;

    const char* s1 = "ATGCTAGCT";
    const char* s2 = "ATGCTAGCT";

    std::vector<int> x1 = convert_dna(s1, mapping);
    std::vector<int> x2 = convert_dna(s2, mapping);

    tree<int> t = construct_tree(x1, 4, 4);

    int matches = compare_to_tree(t, x2, 1);

    std::cout << matches << std::endl;

    // std::vector<int> w;
    // debug(t, 0, w);

    return 0;
}
