#pragma once



template<typename T>
struct vector1D {
    int m_capacity;
    int m_size;
    T *m_data;

    vector1D(int capacity = 0, int size = 0)
        : m_capacity(capacity)
        , m_size(size) {
        m_data = new T[capacity];
    }

    vector1D(const vector1D<T> &oth)
        : vector1D<T>(oth.m_capacity, oth.m_size) {
        std::copy(oth.m_data, oth.m_data + oth.m_size, m_data);
    }

    vector1D(vector1D<T> &&oth) {
        swap(*this, oth);
    }

    friend void swap(vector1D &v1, vector1D &v2) {
        using std::swap;
        swap(v1.m_capacity, v2.m_capacity);
        swap(v1.m_size, v2.m_size);
        swap(v1.m_data, v2.m_data);
    }

    vector1D<T>& operator=(vector1D<T> oth) {
        swap(*this, oth);
        return *this;
    }

    ~vector1D() {
        delete[] m_data;
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

    void push_back(const T& oth) {
        assert(m_size < m_capacity);
        m_data[m_size++] = oth;
    }

    void push_back(T&& oth) {
        assert(m_size < m_capacity);
        m_data[m_size++] = std::move(oth);
    }

    T pop_back() {
        assert(m_size > 0);
        return m_data[--m_size];
    }

    void clear() {
        m_size = 0;
    }

    const T& last() const {
        assert(m_size > 0);
        return m_data[m_size - 1];
    }

    T& last() {
        assert(m_size > 0);
        return m_data[m_size - 1];
    }

    const T& first() const {
        assert(m_size > 0);
        return m_data[0];
    }

    T& first() {
        assert(m_size > 0);
        return m_data[0];
    }


};


template<typename T>
struct vector2D : public vector1D< vector1D<T> > {

    using vector1D< vector1D<T> >::m_data;
    using vector1D< vector1D<T> >::m_capacity;
    using vector1D< vector1D<T> >::m_size;

    vector2D(int rows, int cols)
        : vector1D< vector1D<T> >(rows, rows) {

        for (int i = 0; i < rows; i++) {
            m_data[i] = vector1D<T>(cols);
        }
    }

    ~vector2D() {
    }

};
