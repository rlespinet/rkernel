#pragma once

#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <iostream>

// template<typename T>
// using vector1D = std::vector<T>;

// template<typename T>
// using vector2D = vector1D< vector1D<T> >;

template<typename T>
struct view1D {

    template<typename U>
    friend class vector1D;

    template<typename U>
    friend class vector2D;

    int m_capacity;
    int m_size;
    T *m_data;

    view1D(T* data, int capacity, int size)
        : m_capacity(capacity)
        , m_size(capacity)
        , m_data(data) {
    }

    view1D(T* data = nullptr, int capacity = 0)
        : view1D<T>(data, capacity, capacity) {
    }

    ~view1D() {
    }

    const int capacity() const {
        return m_capacity;
    }

    const int size() const {
        return m_size;
    }

    const T* data() const {
        return m_data;
    }

    T* data() {
        return m_data;
    }

    const T& operator[](int i) const {
        assert(i < m_size);
        return m_data[i];
    }

    T& operator[](int i)  {
        assert(i < m_size);
        return m_data[i];
    }

    const T& operator()(int i) const {
        assert(i < m_size);
        return m_data[i];
    }

    T& operator()(int i)  {
        assert(i < m_size);
        return m_data[i];
    }

    void resize(int new_size) {
        assert(new_size <= capacity);
        this->m_size = new_size;
    }

private:

    view1D(const view1D<T> &oth)
        : m_capacity(oth.m_capacity)
        , m_size(oth.m_size)
        , m_data(oth.m_data) {
    }

    view1D(view1D<T> &&oth) {
        swap(*this, oth);
    }

    view1D<T>& operator=(view1D<T> oth) {
        swap(*this, oth);
        return *this;
    }

    friend void swap(view1D &v1, view1D &v2) {
        using std::swap;
        swap(v1.m_capacity, v2.m_capacity);
        swap(v1.m_size, v2.m_size);
        swap(v1.m_data, v2.m_data);
    }



};


template<typename T>
struct vector1D : public view1D<T> {

    using view1D<T>::m_data;
    using view1D<T>::m_capacity;
    using view1D<T>::m_size;

    vector1D(int capacity = 0)
        : vector1D<T>(capacity, capacity) {
    }

    vector1D(int capacity, int size)
        : view1D<T>(new T[capacity], capacity, size) {
    }

    vector1D(const vector1D<T> &oth)
        : vector1D<T>(oth.m_capacity, oth.m_size) {
        // std::copy(oth.m_data, oth.m_data + oth.m_size, m_data);
        for (int i = 0; i < oth.m_size; i++) {
            m_data[i] = oth.m_data[i];
        }
        assert(0); // Because costly
    }

    vector1D(vector1D<T> &&oth)
        : vector1D<T>() {
        swap(*this, oth);
    }

    ~vector1D() {
        if (m_data != nullptr) {
            delete[] m_data;
        }
    }

    friend void swap(vector1D &v1, vector1D &v2) {
        using std::swap;
        swap(static_cast<view1D<T>&>(v1),static_cast<view1D<T>&>(v2));
    }

    vector1D<T>& operator=(vector1D<T> oth) {
        // TODO(RL) Maybe we don't need to reallocate :/
        swap(*this, oth);
        return *this;
    }

};

template<typename T>
struct vector2D : public vector1D< view1D<T> > {

    using vector1D< view1D<T> >::m_data;
    using vector1D< view1D<T> >::m_capacity;
    using vector1D< view1D<T> >::m_size;

    vector2D(int rows, int cols)
        : vector1D< view1D<T> >(rows) {

        T *data = new T[rows * cols];
        for (int i = 0; i < rows; i++) {
            m_data[i] = view1D<T>(data + i * cols, cols);
        }

    }

    ~vector2D() {
        delete[] (*this)[0].data();
    }

};

template<typename T>
struct matrix {
    T *m_data;
    int m_rows;
    int m_cols;

    matrix()
        : m_data(nullptr)
        , m_rows(0)
        , m_cols(0) {
    }

    matrix(int rows, int cols)
        : m_rows(rows)
        , m_cols(cols) {
        // std::cout << "matrix_alloc" << std::endl;
        m_data = new T[m_rows * m_cols];
    }

    matrix(const matrix& oth) {
        // std::cout << "sq_matrix_copy" << std::endl;
        assert(0);
    }

    matrix(matrix &&oth) {
        // std::cout << "matrix_move" << std::endl;
        swap(*this, oth);
    }

    ~matrix() {
        // std::cout << "matrix_delete" << std::endl;
        if (m_data != nullptr) {
            delete[] m_data;
        }
    }

    matrix& operator=(matrix oth) {
        // std::cout << "sq_matrix_eq" << std::endl;
        swap(*this, oth);
        return *this;
    }

    friend void swap(matrix &m1, matrix &m2) {
        using std::swap;
        // std::cout << "matrix_swap" << std::endl;
        swap(m1.m_rows, m2.m_rows);
        swap(m1.m_cols, m2.m_cols);
        swap(m1.m_data, m2.m_data);
    }

    const T& operator()(int i, int j) const {
        return this->m_data[i * m_cols + j];
    }

    T& operator()(int i, int j) {
        return this->m_data[i * m_cols + j];
    }

    const int rows() const {
        return this->m_rows;
    }

    const int cols() const {
        return this->m_cols;
    }

    const T* data() const {
        return this->m_data;
    }

    T* data() {
        return this->m_data;
    }

    // TODO(RL) Is this a good idea to do this ?
    T* steal_data() {
        T* data = m_data;
        m_data = nullptr;
        this->reset();
        return data;
    }

    void reset() {
        matrix<T> def;
        swap(*this, def);
    }

};

template<typename T>
struct sq_matrix : public matrix<T> {

    sq_matrix(int size)
        : matrix<T>(size, size) {
        // std::cout << "sq_matrix_alloc" << std::endl;
    }

    sq_matrix(const sq_matrix& oth)
        : matrix<T>(oth) {
        // std::cout << "sq_matrix_copy" << std::endl;
        assert(0);
    }

    ~sq_matrix() {
        // std::cout << "sq_matrix_delete" << std::endl;
    }

   friend void swap(sq_matrix &m1, sq_matrix &m2) {
        using std::swap;
        // std::cout << "sq_matrix_swap" << std::endl;
        swap(static_cast<matrix<T>&>(m1), static_cast<matrix<T>&>(m2));
    }

    sq_matrix& operator=(sq_matrix oth) {
        // std::cout << "sq_matrix_eq" << std::endl;
        swap(*this, oth);
        return *this;
    }

    const int size() const {
        return this->m_rows;
    }

    inline void symmetrise_lower() {
        // TODO(RL) This is very cache inefficient
        for (int i = 0; i < size() - 1; i++) {
            for (int j = i + 1; j < size(); j++) {
                this->m_data[j * size() + i] = this->m_data[i * size() + j];
            }
        }
    }

    inline void symmetrise_upper() {
        // TODO(RL) This is very cache inefficient
        for (int i = 0; i < size() - 1; i++) {
            for (int j = i + 1; j < size(); j++) {
                this->m_data[i * size() + j] = this->m_data[j * size() + i];
            }
        }
    }

};
