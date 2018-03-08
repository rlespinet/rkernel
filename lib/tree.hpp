#pragma once

#include "utils.hpp"

template<class T>
class tree {

public:

    typedef int node;

public:

    tree(int children, int max_depth);
    ~tree();

    node parent(node n);

    node child(node n, int i);

    int child_index(node n);

    node first_child(node n);

    bool is_leaf(node n);

    node root();

    T& operator[](node n);
    const T& operator[](node n) const;

    T *data;
    int D;
    int N;
    int size;
};

template<class T>
tree<T>::tree(int children, int max_depth)
    : N(children)
    , D(max_depth) {
    this->size = (ipow(N, max_depth + 1) - 1) / (N - 1);
    this->data = new T[size]();
}


template<class T>
tree<T>::~tree() {
    if (this->data) {
	delete[] this->data;
    }
}

template<class T>
typename tree<T>::node tree<T>::parent(node n) {
    return (n - 1) / N;
}

template<class T>
typename tree<T>::node tree<T>::first_child(node n) {
    return n * N + 1;
}

template<class T>
int tree<T>::child_index(node n) {
    return (n - 1) % N;
}

template<class T>
typename tree<T>::node tree<T>::child(node n, int i) {
    return first_child(n) + i;
}

template<class T>
typename tree<T>::node tree<T>::root() {
    return 0;
}

template<class T>
bool tree<T>::is_leaf(node n) {
    return first_child(n) >= size;
}



template<class T>
T& tree<T>::operator[](node n) {
    return data[n];
}

template<class T>
const T& tree<T>::operator[](node n) const {
    return data[n];
}
