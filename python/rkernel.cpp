#include <Python.h>
#include <iostream>
#include <cstring>
#include <set>
#include <map>
#include <numpy/arrayobject.h>
// #include "rkernel.h"

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
    {NULL, NULL, 0, NULL}
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

inline static int *reindex_array(const int* data, int rows, int cols) {

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
        return NULL;
    }

    for (int i = 0; i < rows * cols; i++) {
        reindexed_data[i] = symbol_map[data[i]];
    }

    return reindexed_data;

}

inline static int *parse_PyArrayString(PyObject *obj, int *prows, int *pcols) {

    PyArray_Descr *des = PyArray_DescrNewFromType(NPY_UNICODE);
    if (des == NULL) {
        return NULL;
    }

    PyObject *str_array = PyArray_FromAny(obj, des, 1, 1,
                                          NPY_ARRAY_C_CONTIGUOUS |
                                          NPY_ARRAY_ALIGNED, NULL);
    if (str_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "Array is of the wrong type");
        // Py_DECREF(des);
        return NULL;
    }

    Py_DECREF(des);

    npy_intp rows = PyArray_DIM(str_array, 0);
    npy_intp cols = PyArray_STRIDE(str_array, 0) / 4;

    int *str_array_data = (int*) PyArray_DATA(str_array);
    if (str_array_data == NULL) {
        PyErr_SetString(PyExc_ValueError, "Bad array format");
        Py_DECREF(str_array);
        return NULL;
    }

    int *int_array_data = reindex_array(str_array_data, rows, cols);
    if (int_array_data == NULL) {
        PyErr_SetString(PyExc_MemoryError, "System ran out of memory");
        Py_DECREF(str_array);
        return NULL;
    }

    Py_DECREF(str_array);

    if (prows != NULL) *prows = rows;
    if (pcols != NULL) *pcols = cols;

    return int_array_data;

}

static PyObject *rkernel_spectrum_bind(PyObject *self, PyObject *args, PyObject* kwargs) {
    PyObject* obj = NULL;

    int k = 3;

    char *keywords[] = {
        "",
        "k",
        NULL
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", keywords, &obj, &k))
        return NULL;

    if (k < 0) {
        PyErr_SetString(PyExc_ValueError, "Substring size k must be positive");
        Py_RETURN_NONE;
    }

    int rows, cols;
    int *array = parse_PyArrayString(obj, &rows, &cols);
    if (array == NULL) {
        Py_RETURN_NONE;
    }

    std::cout << "k = " << k << std::endl;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << array[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }

    delete[] array;
    Py_RETURN_NONE;
}

static PyObject *rkernel_substring_bind(PyObject *self, PyObject *args, PyObject *kwargs) {
    Py_RETURN_NONE;
}


static PyObject *rkernel_mismatch_bind(PyObject *self, PyObject *args, PyObject *kwargs) {
    Py_RETURN_NONE;
}
