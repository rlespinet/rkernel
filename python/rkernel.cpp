#include <Python.h>
#include <iostream>
#include <cstring>

#include "numpy/npy_no_deprecated_api.h"
#include "numpy/arrayobject.h"
#include "rkernel.hpp"

static char module_docstring[] = "Fast kernels for sequences";

static char mismatch_docstring [] = "The mismatch kernel";
static char spectrum_docstring [] = "The spectrum kernel";
static char substring_docstring [] = "The substring kernel";

static PyObject *rkernel_mismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_spectrum_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_substring_bind(PyObject *self, PyObject *args, PyObject *kwargs);

static PyMethodDef rkernel_methods[] = {
    {"spectrum", (PyCFunction) rkernel_spectrum_bind, METH_VARARGS | METH_KEYWORDS, spectrum_docstring},
    {"substring", (PyCFunction) rkernel_substring_bind, METH_VARARGS  | METH_KEYWORDS, substring_docstring},
    {"mismatch", (PyCFunction) rkernel_mismatch_bind, METH_VARARGS  | METH_KEYWORDS, mismatch_docstring},
    {nullptr, nullptr, 0, nullptr}
};

static struct PyModuleDef rkernel_module = {
    PyModuleDef_HEAD_INIT,
    "rkernel",
    module_docstring,
    -1,
    rkernel_methods
};

PyMODINIT_FUNC PyInit_rkernel(void) {
    import_array();

    return PyModule_Create(&rkernel_module);
}

template<typename letter>
inline static sequence_array<letter> parse_PyArrayString(PyObject *obj) {

    PyArrayObject *str_array = (PyArrayObject*) PyArray_FROMANY(obj, NPY_UNICODE, 1, 1,
                                                                NPY_ARRAY_C_CONTIGUOUS |
                                                                NPY_ARRAY_ALIGNED);
    if (str_array == nullptr) {
        PyErr_SetString(PyExc_TypeError, "Array is of the wrong type");
        throw "TODO";
    }

    npy_intp rows = PyArray_DIM(str_array, 0);
    npy_intp cols = PyArray_STRIDE(str_array, 0) / 4;

    int *str_array_data = (int*) PyArray_DATA(str_array);
    if (str_array_data == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Bad array format");
        Py_DECREF(str_array);
        throw "TODO";
    }

    sequence_array<letter> seq_array = to_sequence_array<letter, int>(str_array_data, rows, cols);

    Py_DECREF(str_array);

    return seq_array;
}

static PyObject *rkernel_spectrum_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;

    char *keywords[] = {"", "k", nullptr};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", keywords, &obj, &k))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    sequence_array<int> seq_array = parse_PyArrayString<int>(obj);

    sq_matrix<float> kernel = spectrum<int, float>(seq_array.sequences,
                                                   seq_array.sequences_len,
                                                   seq_array.alphabet_size, k);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}

static PyObject *rkernel_substring_bind(PyObject *self, PyObject *args, PyObject *kwargs) {
    Py_RETURN_NONE;
}

static PyObject *rkernel_mismatch_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    Py_RETURN_NONE;
}

// static PyObject *rkernel_mismatch_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
//     PyObject* obj = nullptr;

//     int k = 3;
//     int m = 1;

//     char *keywords[] = {
//         "",
//         "k",
//         "m",
//         nullptr
//     };

//     if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ii", keywords, &obj, &k, &m))
//         return nullptr;

//     if (k < 0) {
//         PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
//         Py_RETURN_NONE;
//     }

//     if (m < 0) {
//         PyErr_SetString(PyExc_ValueError, "Allowed mismatch must be positive");
//         Py_RETURN_NONE;
//     }

//     sequence_array<int> seq_array = parse_PyArrayString<int>(obj);

//     sq_matrix<float> kernel = mismatch<int, float>(seq_array.sequences,
//                                                    seq_array.sequences_len,
//                                                    seq_array.alphabet_size, k, m);

//     npy_intp dims[] = {kernel.rows(), kernel.cols()};
//     PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
//     if (kernel_matrix == nullptr) {
//         PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
//         Py_RETURN_NONE;
//     }

//     return kernel_matrix;

// }
