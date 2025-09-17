#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "symnmf.h"

static PyObject *fit(PyObject *self, PyObject *args);

static PyObject *sym(PyObject *self, PyObject *args)
{
    PyObject *data_py;
    const char *goal;
    if (!PyArg_ParseTuple(args, "O", &data_py))
    {
        Py_RETURN_NONE;
    }
    goal = "sym";
    return fit(self, Py_BuildValue("sO", goal, data_py));
}

static PyObject *ddg(PyObject *self, PyObject *args)
{
    PyObject *data_py;
    const char *goal;
    if (!PyArg_ParseTuple(args, "O", &data_py))
    {
        Py_RETURN_NONE;
    }
    goal = "ddg";
    return fit(self, Py_BuildValue("sO", goal, data_py));
}

static PyObject *norm(PyObject *self, PyObject *args)
{
    PyObject *data_py;
    const char *goal;
    if (!PyArg_ParseTuple(args, "O", &data_py))
    {
        Py_RETURN_NONE;
    }
    goal = "norm";
    return fit(self, Py_BuildValue("sO", goal, data_py));
}

// Serve the functions sym, ddg and norm to the python file
static PyObject *fit(PyObject *self, PyObject *args)
{
    PyObject *data_py;
    const char *goal;
    Py_ssize_t n, D;
    PyObject *row_obj;
    int i, j, N, d;
    double **sol;
    double ***sol2;
    double *arr;
    double num;
    if (!PyArg_ParseTuple(args, "sO", &goal, &data_py))
    {
        Py_RETURN_NONE;
    }
    n = PyList_Size(data_py);
    N = (int)n;
    row_obj = PyList_GetItem(data_py, 0);
    D = PyList_Size(row_obj);
    d = (int)D;
    double **data = calloc(N, sizeof(double *));
    if (data == NULL)
    {
        Py_RETURN_NONE;
    }
    for (i = 0; i < N; i++)
    {
        data[i] = calloc(d, sizeof(double));
        if (data[i] == NULL)
        {
            Py_RETURN_NONE;
        }
        for (j = 0; j < d; j++)
        {
            PyObject *row_obj = PyList_GetItem(data_py, i);
            PyObject *item = PyList_GetItem(row_obj, j);
            num = PyFloat_AsDouble(item);
            data[i][j] = num;
        }
    }
    if (strcmp(goal, "sym") == 0)
    {
        sol = symC(data, N, d);
        if (sol == NULL) // if sym failed
        {
            Py_RETURN_NONE;
        }
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        sol2 = ddgC(data, N, d);
        if (sol2 == NULL) // if ddg failed
        {
            Py_RETURN_NONE;
        }
        sol = sol2[1];
    }
    else
    {
        sol = normC(data, N, d);
        if (sol == NULL) // if norm failed
        {
            Py_RETURN_NONE;
        }
    }
    arr = calloc(N * N, sizeof(double));
    if (arr == NULL)
    {
        Py_RETURN_NONE;
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            arr[(i * N) + j] = sol[i][j];
        }
    }
    PyObject *result = PyTuple_New(N * N);
    for (Py_ssize_t i = 0; i < N * N; i++)
    {
        PyTuple_SetItem(result, i, PyFloat_FromDouble(arr[i]));
    }
    return result;
}

// Serve the full symNMF algorithm to the python file
static PyObject *symnmf(PyObject *self, PyObject *args)
{
    PyObject *H_py, *W_py;
    int N, k, i, j;
    double **sol;
    double *arr;
    double num;
    if (!PyArg_ParseTuple(args, "iiOO", &N, &k, &H_py, &W_py))
    {
        Py_RETURN_NONE;
    }
    double **H = calloc(N, sizeof(double *));
    if (H == NULL)
    {
        Py_RETURN_NONE;
    }
    for (i = 0; i < N; i++)
    {
        H[i] = calloc(k, sizeof(double));
        if (H[i] == NULL)
        {
            Py_RETURN_NONE;
        }
        for (j = 0; j < k; j++)
        {
            PyObject *item = PyList_GetItem(H_py, (k * i) + j);
            num = PyFloat_AsDouble(item);
            H[i][j] = num;
        }
    }
    double **W = calloc(N, sizeof(double *));
    if (W == NULL)
    {
        Py_RETURN_NONE;
    }
    for (i = 0; i < N; i++)
    {
        W[i] = calloc(N, sizeof(double));
        if (W[i] == NULL)
        {
            Py_RETURN_NONE;
        }
        for (j = 0; j < N; j++)
        {
            PyObject *item = PyList_GetItem(W_py, (N * i) + j);
            num = PyFloat_AsDouble(item);
            W[i][j] = num;
        }
    }
    sol = symnmfC(H, W, N, k);
    arr = calloc(N * k, sizeof(double));
    if (arr == NULL)
    {
        Py_RETURN_NONE;
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < k; j++)
        {
            arr[(i * k) + j] = sol[i][j];
        }
    }
    PyObject *result = PyTuple_New(N * k);
    for (Py_ssize_t i = 0; i < N * k; i++)
    {
        PyTuple_SetItem(result, i, PyFloat_FromDouble(arr[i]));
    }
    return result;
}

static PyMethodDef symnmfMethods[] = {
    {"fit",
     (PyCFunction)fit,
     METH_VARARGS,
     PyDoc_STR("The function expects *data and returns matrix according to goal")},
    {"sym",
     (PyCFunction)sym,
     METH_VARARGS,
     PyDoc_STR("The function expects *data and returns matrix to fit according to the goal received")},
    {"ddg",
     (PyCFunction)ddg,
     METH_VARARGS,
     PyDoc_STR("The function expects *data and returns matrix to fit according to the goal received")},
    {"norm",
     (PyCFunction)norm,
     METH_VARARGS,
     PyDoc_STR("The function expects *data and returns matrix to fit according to the goal received")},
    {"symnmf",
     (PyCFunction)symnmf,
     METH_VARARGS,
     PyDoc_STR("The function expects int N, int k, *H, *W and returns the clustering matrix")},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    NULL,
    -1,
    symnmfMethods};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m)
    {
        Py_RETURN_NONE;
    }
    return m;
}