#pragma once

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
