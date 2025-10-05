#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "sasa_core.h"

static PyObject *compute_sasa_wrapper(PyObject *self, PyObject *args, PyObject *kwargs) {
    PyObject *coords_obj = NULL;
    PyObject *radii_obj = NULL;
    PyArrayObject *coords_array = NULL;
    PyArrayObject *radii_array = NULL;
    float probe_radius = 1.4f;
    int n_samples = 500;
    unsigned int seed = VMD_SEED;

    static char *kwlist[] = {"coords", "radii", "probe_radius", "n_samples", "seed", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|fiI", kwlist,
                                     &coords_obj, &radii_obj,
                                     &probe_radius, &n_samples, &seed)) {
        return NULL;
    }

    // Validate parameters
    if (n_samples <= 0) {
        PyErr_SetString(PyExc_ValueError, "n_samples must be positive");
        return NULL;
    }

    // Convert input objects to numpy arrays (allow any numeric type)
    coords_array = (PyArrayObject*)PyArray_FromAny(coords_obj, NULL,
                                                   2, 2, NPY_ARRAY_CARRAY_RO, NULL);
    if (coords_array == NULL) {
        return NULL;
    }

    radii_array = (PyArrayObject*)PyArray_FromAny(radii_obj, NULL,
                                                  1, 1, NPY_ARRAY_CARRAY_RO, NULL);
    if (radii_array == NULL) {
        Py_DECREF(coords_array);
        return NULL;
    }

    // Validate input arrays
    if (PyArray_NDIM(coords_array) != 2 || PyArray_DIM(coords_array, 1) != 3) {
        PyErr_SetString(PyExc_ValueError, "coords must be a 2D array with shape (n_atoms, 3)");
        Py_DECREF(coords_array);
        Py_DECREF(radii_array);
        return NULL;
    }

    if (PyArray_NDIM(radii_array) != 1) {
        PyErr_SetString(PyExc_ValueError, "radii must be a 1D array");
        Py_DECREF(coords_array);
        Py_DECREF(radii_array);
        return NULL;
    }

    int n_atoms = (int)PyArray_DIM(coords_array, 0);
    if (PyArray_DIM(radii_array, 0) != n_atoms) {
        PyErr_SetString(PyExc_ValueError, "coords and radii arrays must have same number of atoms");
        Py_DECREF(coords_array);
        Py_DECREF(radii_array);
        return NULL;
    }

    // Convert to float arrays and ensure contiguous
    PyArrayObject *coords_float = (PyArrayObject*)PyArray_Cast(coords_array, NPY_FLOAT32);
    PyArrayObject *radii_float = (PyArrayObject*)PyArray_Cast(radii_array, NPY_FLOAT32);

    Py_DECREF(coords_array);
    Py_DECREF(radii_array);

    if (!coords_float || !radii_float) {
        Py_XDECREF(coords_float);
        Py_XDECREF(radii_float);
        return NULL;
    }

    // Ensure arrays are contiguous
    PyArrayObject *coords_contig = (PyArrayObject*)PyArray_GETCONTIGUOUS(coords_float);
    PyArrayObject *radii_contig = (PyArrayObject*)PyArray_GETCONTIGUOUS(radii_float);

    Py_DECREF(coords_float);
    Py_DECREF(radii_float);

    if (!coords_contig || !radii_contig) {
        Py_XDECREF(coords_contig);
        Py_XDECREF(radii_contig);
        return NULL;
    }

    // Get data pointers
    float *coords_data = (float*)PyArray_DATA(coords_contig);
    float *radii_data = (float*)PyArray_DATA(radii_contig);

    // Validate for NaN and infinite values
    for (int i = 0; i < n_atoms * 3; i++) {
        if (!isfinite(coords_data[i])) {
            PyErr_SetString(PyExc_ValueError, "coords contains NaN or infinite values");
            Py_DECREF(coords_contig);
            Py_DECREF(radii_contig);
            return NULL;
        }
    }

    for (int i = 0; i < n_atoms; i++) {
        if (!isfinite(radii_data[i])) {
            PyErr_SetString(PyExc_ValueError, "radii contains NaN or infinite values");
            Py_DECREF(coords_contig);
            Py_DECREF(radii_contig);
            return NULL;
        }
    }

    // Call C implementation
    float total_sasa;
    PointList surface_points = {0};

    int result = compute_sasa(coords_data, radii_data, n_atoms,
                             probe_radius, n_samples, seed,
                             &total_sasa, &surface_points);

    Py_DECREF(coords_contig);
    Py_DECREF(radii_contig);

    if (result != 0) {
        free_point_list(&surface_points);
        PyErr_SetString(PyExc_RuntimeError, "SASA computation failed");
        return NULL;
    }

    // Create numpy array for surface points
    npy_intp dims[2] = {surface_points.count, 3};
    PyArrayObject *points_array = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_FLOAT32);
    if (!points_array) {
        free_point_list(&surface_points);
        return NULL;
    }

    // Copy surface points to numpy array
    float *points_data = (float*)PyArray_DATA(points_array);
    for (int i = 0; i < surface_points.count; i++) {
        points_data[3*i] = surface_points.points[i].x;
        points_data[3*i+1] = surface_points.points[i].y;
        points_data[3*i+2] = surface_points.points[i].z;
    }

    free_point_list(&surface_points);

    // Return tuple (total_sasa, surface_points)
    PyObject *sasa_float = PyFloat_FromDouble((double)total_sasa);
    if (!sasa_float) {
        Py_DECREF(points_array);
        return NULL;
    }

    PyObject *result_tuple = PyTuple_New(2);
    if (!result_tuple) {
        Py_DECREF(sasa_float);
        Py_DECREF(points_array);
        return NULL;
    }

    PyTuple_SET_ITEM(result_tuple, 0, sasa_float);
    PyTuple_SET_ITEM(result_tuple, 1, (PyObject*)points_array);

    return result_tuple;
}

static PyMethodDef sasa_methods[] = {
    {"compute_sasa", (PyCFunction)compute_sasa_wrapper, METH_VARARGS | METH_KEYWORDS,
     "Compute SASA using Monte Carlo sampling.\n\n"
     "Parameters:\n"
     "  coords: numpy array of shape (n_atoms, 3) with atom coordinates\n"
     "  radii: numpy array of shape (n_atoms,) with atom radii\n"
     "  probe_radius: float, probe radius (default 1.4)\n"
     "  n_samples: int, number of Monte Carlo samples per atom (default 500)\n"
     "  seed: int, random seed for reproducibility (default 38572111)\n\n"
     "Returns:\n"
     "  tuple: (total_sasa, surface_points)\n"
     "    total_sasa: float, total solvent accessible surface area\n"
     "    surface_points: numpy array of shape (n_points, 3) with surface point coordinates"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef sasa_module = {
    PyModuleDef_HEAD_INIT,
    "sasa_ext",
    "Fast SASA computation using Monte Carlo sampling",
    -1,
    sasa_methods
};

PyMODINIT_FUNC PyInit_sasa_ext(void) {
    import_array();  // Initialize numpy
    return PyModule_Create(&sasa_module);
}