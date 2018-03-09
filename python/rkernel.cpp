#include <Python.h>
#include <iostream>
#include <cstring>

#include "numpy/npy_no_deprecated_api.h"
#include "numpy/arrayobject.h"
#include "rkernel.hpp"

static char module_docstring[] = "Fast kernels for sequence";

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

inline static seq_array_t parse_PyArrayString(PyObject *obj) {

    PyArray_Descr *des = PyArray_DescrNewFromType(NPY_UNICODE);
    if (des == nullptr) {
        return invalid_seq_array;
    }

    PyArrayObject *str_array = (PyArrayObject*) PyArray_FromAny(obj, des, 1, 1,
                                          NPY_ARRAY_C_CONTIGUOUS |
                                          NPY_ARRAY_ALIGNED, nullptr);
    if (str_array == nullptr) {
        PyErr_SetString(PyExc_TypeError, "Array is of the wrong type");
        // Py_DECREF(des);
        return invalid_seq_array;
    }

    // Py_DECREF(des);

    npy_intp rows = PyArray_DIM(str_array, 0);
    npy_intp cols = PyArray_STRIDE(str_array, 0) / 4;

    int *str_array_data = (int*) PyArray_DATA(str_array);
    if (str_array_data == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Bad array format");
        Py_DECREF(str_array);
        return invalid_seq_array;
    }

    seq_array_t seq_array = format_sequence_array(str_array_data, rows, cols);
    if (is_invalid(seq_array)) {
        PyErr_SetString(PyExc_ValueError, "Failed to create the sequence array");
        Py_DECREF(str_array);
        return invalid_seq_array;
    }

    Py_DECREF(str_array);

    return seq_array;
}

static PyObject *rkernel_spectrum_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = nullptr;

    int k = 3;

    char *keywords[] = {
        "",
        "k",
        nullptr
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", keywords, &obj, &k))
        return nullptr;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    seq_array_t seq_array = parse_PyArrayString(obj);
    if (is_invalid(seq_array)) {
        Py_RETURN_NONE;
    }

    kernel_t kernel = spectrum(seq_array, k);
    if (is_invalid(kernel)) {
        PyErr_SetString(PyExc_ValueError, "Failed to compute spectrum kernel");
        Py_RETURN_NONE;
    }

    npy_intp dims[] = {kernel.size, kernel.size};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.data);
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;

}

static PyObject *rkernel_substring_bind(PyObject *self, PyObject *args, PyObject *kwargs) {
    Py_RETURN_NONE;
}


static PyObject *rkernel_mismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs) {
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

    seq_array_t seq_array = parse_PyArrayString(obj);
    if (is_invalid(seq_array)) {
        Py_RETURN_NONE;
    }

    kernel_t kernel = mismatch(seq_array, k, m);
    if (is_invalid(kernel)) {
        PyErr_SetString(PyExc_ValueError, "Failed to compute spectrum kernel");
        Py_RETURN_NONE;
    }

    npy_intp dims[] = {kernel.size, kernel.size};
    PyObject* kernel_matrix = PyArray_SimpleNewFromData(2, dims, NPY_FLOAT, kernel.data);
    if (kernel_matrix == nullptr) {
        PyErr_SetString(PyExc_ValueError, "Failed to construct the final matrix");
        Py_RETURN_NONE;
    }

    return kernel_matrix;


    Py_RETURN_NONE;
}
