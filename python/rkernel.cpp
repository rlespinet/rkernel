#include <Python.h>
#include <iostream>
#include <cstring>

#include "numpy/npy_no_deprecated_api.h"
#include "numpy/arrayobject.h"
#include "rkernel.hpp"

static char module_docstring[] = "Fast kernels for sequences";

static char mismatch_docstring [] = "The mismatch kernel";
static char wmismatch_docstring [] = "The weighted mismatch kernel";
static char dimismatch_docstring [] = "The dimismatch kernel";
static char rmismatch_docstring [] = "The rmismatch kernel";
static char spectrum_docstring [] = "The spectrum kernel";
static char substring_docstring [] = "The substring kernel";
static char wgappy_docstring [] = "The weighted gappy kernel";
static char wildcard_docstring [] = "The wildcard kernel";
static char smith_waterman_docstring [] = "Smith waterman pairwise score";
static char needleman_wunsch_docstring [] = "Smith waterman pairwise score";

static PyObject *rkernel_mismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_wmismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_dimismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_rmismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_spectrum_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_substring_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_wgappy_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_wildcard_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_smith_waterman_bind(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *rkernel_needleman_wunsch_bind(PyObject *self, PyObject *args, PyObject *kwargs);

static PyMethodDef rkernel_methods[] = {
    {"spectrum", (PyCFunction) rkernel_spectrum_bind, METH_VARARGS | METH_KEYWORDS, spectrum_docstring},
    {"substring", (PyCFunction) rkernel_substring_bind, METH_VARARGS  | METH_KEYWORDS, substring_docstring},
    {"mismatch", (PyCFunction) rkernel_mismatch_bind, METH_VARARGS  | METH_KEYWORDS, mismatch_docstring},
    {"wmismatch", (PyCFunction) rkernel_wmismatch_bind, METH_VARARGS  | METH_KEYWORDS, wmismatch_docstring},
    {"dimismatch", (PyCFunction) rkernel_dimismatch_bind, METH_VARARGS  | METH_KEYWORDS, dimismatch_docstring},
    {"rmismatch", (PyCFunction) rkernel_rmismatch_bind, METH_VARARGS  | METH_KEYWORDS, rmismatch_docstring},
    {"wgappy", (PyCFunction) rkernel_wgappy_bind, METH_VARARGS  | METH_KEYWORDS, wgappy_docstring},
    {"wildcard", (PyCFunction) rkernel_wildcard_bind, METH_VARARGS  | METH_KEYWORDS, wildcard_docstring},
    {"needleman_wunsch", (PyCFunction) rkernel_needleman_wunsch_bind, METH_VARARGS  | METH_KEYWORDS, needleman_wunsch_docstring},
    {"smith_waterman", (PyCFunction) rkernel_smith_waterman_bind, METH_VARARGS  | METH_KEYWORDS, smith_waterman_docstring},
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

inline static sequence_array parse_PyArrayString(PyObject *obj) {

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

    sequence_array seq_array = to_sequence_array<int>(str_array_data, rows, cols);

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

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = spectrum<float>(seq_array.sequences,
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
    PyObject* obj = nullptr;

    int k = 3;
    int m = 1;

    char *keywords[] = {
        "",
        "k",
        "m",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ii", keywords, &obj, &k, &m))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "Allowed mismatch must be positive");
        Py_RETURN_NONE;
    }

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = mismatch<float>(seq_array.sequences,
					      seq_array.sequences_len,
					      seq_array.alphabet_size, k, m);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}

static PyObject *rkernel_wmismatch_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;
    int m = 1;
    float w = 1.0f;

    char *keywords[] = {
        "",
        "k",
        "m",
        "w",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iif", keywords, &obj, &k, &m, &w))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "Allowed mismatch must be positive");
        Py_RETURN_NONE;
    }

    if (w < 0.0f) {
        PyErr_SetString(PyExc_ValueError, "Weight must be positive");
        Py_RETURN_NONE;
    }

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = wmismatch<float>(seq_array.sequences,
                                               seq_array.sequences_len,
                                               seq_array.alphabet_size, k, m, w);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}


static PyObject *rkernel_dimismatch_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;
    int m = 1;

    char *keywords[] = {
        "",
        "k",
        "m",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ii", keywords, &obj, &k, &m))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "Allowed mismatch must be positive");
        Py_RETURN_NONE;
    }

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = dimismatch<float>(seq_array.sequences,
                                                seq_array.sequences_len,
                                                seq_array.alphabet_size, k, m);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;
}

static PyObject *rkernel_rmismatch_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;
    int m = 1;
    int r = 2;

    char *keywords[] = {
        "",
        "k",
        "m",
        "r",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iii", keywords, &obj, &k, &m, &r))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "Allowed mismatch must be positive");
        Py_RETURN_NONE;
    }

    if (r <= 0) {
        PyErr_SetString(PyExc_ValueError, "Allowed order must be greater than 1");
        Py_RETURN_NONE;
    }

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = rmismatch<float>(seq_array.sequences,
                                               seq_array.sequences_len,
                                               seq_array.alphabet_size, k, m, r);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;
}

static PyObject *rkernel_wgappy_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;
    int g = 5;
    float w = 1.0f;

    char *keywords[] = {
        "",
        "k",
        "g",
        "w",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iif", keywords, &obj, &k, &g, &w))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    if (g < k) {
        PyErr_SetString(PyExc_ValueError, "Substring size g must be higher than k");
        Py_RETURN_NONE;
    }

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = wgappy<float>(seq_array.sequences,
                                            seq_array.sequences_len,
                                            seq_array.alphabet_size,
                                            k, g, w);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}


static PyObject *rkernel_wildcard_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;
    int m = 1;
    float w = 1.0f;

    char *keywords[] = {
        "",
        "k",
        "m",
        "w",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iif", keywords, &obj, &k, &m, &w))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    if (m < 0) {
        PyErr_SetString(PyExc_ValueError, "Maximum number of wildcard must be positive");
        Py_RETURN_NONE;
    }

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = wildcard<float>(seq_array.sequences,
                                              seq_array.sequences_len,
                                              seq_array.alphabet_size,
                                              k, m, w);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}


static PyObject *rkernel_smith_waterman_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    float gap_open = 10.0f;
    float gap_extension = 0.5f;
    float similarity = 5.0f;
    float substitution = 4.0f;

    char *keywords[] = {
        "",
        "gap_open",
        "gap_extension",
        "similarity",
        "substitution",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ffff", keywords, &obj, &gap_open, &gap_extension, &similarity, &substitution))
        return nullptr;


    // TODO(RL) Verify arguments

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = smith_waterman_affine_gap<float>(seq_array.sequences,
                                                               seq_array.sequences_len,
                                                               seq_array.alphabet_size,
                                                               similarity, substitution,
                                                               gap_open, gap_extension);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}

static PyObject *rkernel_needleman_wunsch_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    float gap_open = 10.0f;
    float gap_extension = 0.5f;
    float similarity = 5.0f;
    float substitution = 4.0f;

    char *keywords[] = {
        "",
        "gap_open",
        "gap_extension",
        "similarity",
        "substitution",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|ffff", keywords, &obj, &gap_open, &gap_extension, &similarity, &substitution))
        return nullptr;


    // TODO(RL) Verify arguments

    sequence_array seq_array = parse_PyArrayString(obj);

    sq_matrix<float> kernel = needleman_wunsch_affine_gap<float>(seq_array.sequences,
                                                                 seq_array.sequences_len,
                                                                 seq_array.alphabet_size,
                                                                 similarity, substitution,
                                                                 gap_open, gap_extension);

    npy_intp dims[] = {kernel.rows(), kernel.cols()};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.steal_data());
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}
