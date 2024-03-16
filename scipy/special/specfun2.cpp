#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

#include "special/specfun.h"

// This is a simple generic inline function that generates the strided loop needed for a ufunc given a function pointer.
// It should go in a common header.
// This can be done for an arbitrary number of arguments.
template <typename Arg0, typename Res, Res (*F)(Arg0)>
void SpecFun_UFuncLoop(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
    for (npy_intp i = 0; i < dimensions[0]; ++i) {
        *reinterpret_cast<Res *>(args[1]) = F(*reinterpret_cast<Arg0 *>(args[0]));

        args[0] += steps[0];
        args[1] += steps[1];
    }
}

// Just initializes everything needed, can also go in a common header
bool SpecFun_Initialize() {
    Py_Initialize();

    import_array();
    if (PyErr_Occurred()) {
        // import array failed
        return false;
    }

    import_umath();
    if (PyErr_Occurred()) {
        // import umath failed
        return false;
    }

    return true;
}

#if (PY_VERSION_HEX < 0x030a00f0)
// this is in Python >=3.10
int PyModule_AddObjectRef(PyObject *module, const char *name, PyObject *value) {
    Py_INCREF(value);
    return PyModule_AddObject(module, name, value);
}
#endif

// These two definitions will be needed for each ufunc
static PyUFuncGenericFunction expi_funcs[2] = {
    SpecFun_UFuncLoop<double, double, special::expi>,
    SpecFun_UFuncLoop<std::complex<double>, std::complex<double>, special::cexpi>};
static char expi_types[2][2] = {{NPY_FLOAT64, NPY_FLOAT64}, {NPY_COMPLEX128, NPY_COMPLEX128}};

static PyModuleDef specfun2_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "specfun2",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_specfun2() {
    SpecFun_Initialize();

    PyObject *specfun2 = PyModule_Create(&specfun2_def);
    if (specfun2 == nullptr) {
        return nullptr;
    }

    // These two lines will be needed for each ufunc
    PyObject *expi = PyUFunc_FromFuncAndData(expi_funcs, nullptr, reinterpret_cast<char *>(expi_types), 2, 1, 1,
                                             PyUFunc_None, "expi", nullptr, 0);
    PyModule_AddObjectRef(specfun2, "expi", expi);

    return specfun2;
}