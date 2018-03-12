#pragma once

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
